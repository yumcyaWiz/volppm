#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include <omp.h>

#include <optional>

#include "camera.h"
#include "core.h"
#include "image.h"
#include "photon_map.h"
#include "scene.h"

class Integrator {
 protected:
  const std::shared_ptr<Camera> camera;

 public:
  Integrator(const std::shared_ptr<Camera>& camera) : camera(camera) {}

  // render scene
  virtual void render(const Scene& scene, Sampler& sampler, Image& image) = 0;

  // compute cosine term
  // NOTE: need to account for the asymmetry of BSDF when photon tracing
  // https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#x3-Non-symmetryDuetoShadingNormals
  // Veach, Eric. Robust Monte Carlo methods for light transport simulation.
  // Stanford University, 1998. Section 5.3
  static float cosTerm(const Vec3f& wo, const Vec3f& wi,
                       const SurfaceInfo& surfaceInfo,
                       const TransportDirection& transport_dir) {
    const float wi_ns = dot(wi, surfaceInfo.shadingNormal);
    const float wi_ng = dot(wi, surfaceInfo.geometricNormal);
    const float wo_ns = dot(wo, surfaceInfo.shadingNormal);
    const float wo_ng = dot(wo, surfaceInfo.geometricNormal);

    // prevent light leaks
    if (wi_ng * wi_ns <= 0 || wo_ng * wo_ns <= 0) {
      return 0;
    }

    if (transport_dir == TransportDirection::FROM_CAMERA) {
      return std::abs(wi_ns);
    } else if (transport_dir == TransportDirection::FROM_LIGHT) {
      return std::abs(wo_ns) * std::abs(wi_ng) / std::abs(wo_ng);
    } else {
      spdlog::error("[Integrator] invalid transport direction");
      std::exit(EXIT_FAILURE);
    }
  }

  static bool isTransmitted(const Vec3f& wo, const Vec3f& wi, const Vec3f& n) {
    return (dot(wo, n) < 0) != (dot(wi, n) < 0);
  }

  static bool isEntered(const Vec3f& wi, const Vec3f& n) {
    return dot(wi, n) < 0;
  }

  // push or pop medium
  static void updateMedium(Ray& ray, const Vec3f& wi,
                           const IntersectInfo& info) {
    if (isTransmitted(-ray.direction, wi, info.surfaceInfo.shadingNormal)) {
      if (isEntered(wi, info.surfaceInfo.shadingNormal)) {
        if (info.hitPrimitive->hasMedium()) {
          ray.pushMedium(info.hitPrimitive->getMedium());
        }
      } else {
        ray.popMedium();
      }
    }
  }
};

// abstraction of path based integrator
class PathIntegrator : public Integrator {
 private:
  // number of samples in each pixel
  const uint32_t n_samples;

 public:
  // compute radiance coming from the given ray
  virtual Vec3f integrate(const Ray& ray, const Scene& scene,
                          Sampler& sampler) const = 0;

  PathIntegrator(const std::shared_ptr<Camera>& camera, uint32_t n_samples)
      : Integrator(camera), n_samples(n_samples) {}

  void render(const Scene& scene, Sampler& sampler,
              Image& image) override final {
    const uint32_t width = image.getWidth();
    const uint32_t height = image.getHeight();

    spdlog::info("[PathIntegrator] rendering...");
#pragma omp parallel for collapse(2) schedule(dynamic, 1)
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        // init sampler for each pixel
        const std::unique_ptr<Sampler> sampler_per_pixel = sampler.clone();
        sampler_per_pixel->setSeed((sampler.getSeed() + 1) * (j + width * i));

        // warmup sampler
        for (uint32_t k = 0; k < 100; ++k) {
          sampler_per_pixel->getNext1D();
        }

        // iteration
        for (uint32_t k = 0; k < n_samples; ++k) {
          // SSAA
          const float u =
              -(2.0f * (j + sampler_per_pixel->getNext1D()) - width) / height;
          const float v =
              (2.0f * (i + sampler_per_pixel->getNext1D()) - height) / height;

          Ray ray;
          float pdf;
          if (camera->sampleRay(Vec2f(u, v), *sampler_per_pixel, ray, pdf)) {
            // compute incoming radiance
            const Vec3f radiance =
                integrate(ray, scene, *sampler_per_pixel) / pdf;

            // invalid radiance check
            if (std::isnan(radiance[0]) || std::isnan(radiance[1]) ||
                std::isnan(radiance[2])) {
              spdlog::error("[PathIntegrator] radiance is NaN");
              continue;
            } else if (std::isinf(radiance[0]) || std::isinf(radiance[1]) ||
                       std::isinf(radiance[2])) {
              spdlog::error("[PathIntegrator] radiance is inf");
              continue;
            } else if (radiance[0] < 0 || radiance[1] < 0 || radiance[2] < 0) {
              spdlog::error("[PathIntegrator] radiance is minus");
              continue;
            }

            image.addPixel(i, j, radiance);
          } else {
            image.setPixel(i, j, Vec3f(0));
          }
        }
      }
    }
    spdlog::info("[PathIntegrator] done");

    // take average
    image /= Vec3f(n_samples);
  }
};

// implementation of volumetric unidirectional path tracing
class PathTracing : public PathIntegrator {
 private:
  const uint32_t maxDepth;

 public:
  PathTracing(const std::shared_ptr<Camera>& camera, uint32_t n_samples,
              uint32_t maxDepth = 100)
      : PathIntegrator(camera, n_samples), maxDepth(maxDepth) {}

  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const override {
    Vec3f radiance(0);
    Ray ray = ray_in;
    ray.throughput = Vec3f(1, 1, 1);

    uint32_t depth = 0;
    while (depth < maxDepth) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (depth > 0) {
          const float russian_roulette_prob = std::min(
              (ray.throughput[0] + ray.throughput[1] + ray.throughput[2]) /
                  3.0f,
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          ray.throughput /= Vec3f(russian_roulette_prob);
        }

        // sample medium
        bool is_scattered = false;
        if (ray.hasMedium()) {
          const Medium* medium = ray.getCurrentMedium();

          Vec3f pos;
          Vec3f dir;
          Vec3f throughput_medium;
          is_scattered = medium->sampleMedium(ray, info.t, sampler, pos, dir,
                                              throughput_medium);

          // advance ray
          ray.origin = pos;
          ray.direction = dir;

          // update throughput
          ray.throughput *= throughput_medium;
        }

        bool is_reflected_or_refracted = false;
        if (!is_scattered) {
          // Le
          if (info.hitPrimitive->hasAreaLight()) {
            radiance += ray.throughput *
                        info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
            break;
          }

          // sample direction by BxDF
          Vec3f dir = ray.direction;
          if (info.hitPrimitive->hasSurface()) {
            float pdf_dir;
            const Vec3f f = info.hitPrimitive->sampleBxDF(
                -ray.direction, info.surfaceInfo,
                TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

            // update throughput
            ray.throughput *= f *
                              cosTerm(-ray.direction, dir, info.surfaceInfo,
                                      TransportDirection::FROM_CAMERA) /
                              pdf_dir;

            is_reflected_or_refracted = true;
          }

          // update ray's medium
          updateMedium(ray, dir, info);

          // update ray
          ray.origin = info.surfaceInfo.position;
          ray.direction = dir;
        }

        // update depth
        if (is_scattered || is_reflected_or_refracted) {
          depth++;
        }
      } else {
        // ray goes out to the sky
        break;
      }
    }

    return radiance;
  }
};

// this implementation is based on modified version of original PPM
//  Knaus, Claude, and Matthias Zwicker.
// "Progressive photon mapping: A probabilistic approach." ACM Transactions on
// Graphics (TOG) 30.3 (2011): 1-13.
class PPM : public Integrator {
 private:
  // number of iterations
  const uint32_t nIterations;
  // number of photons in each iteration
  const uint32_t nPhotons;
  // parameter for radius reduction, see the paper
  const float alpha;
  // maximum tracing depth
  const uint32_t maxDepth;

  // number of emitted photons
  uint32_t nEmittedPhotons;
  // global photon search radius for surface radiance estimation
  float globalRadius;
  // global photon search radius for volume radiance estimation
  float globalRadiusVolume;

  PhotonMap photonMap;
  PhotonMap volumePhotonMap;

  // compute reflected radiance with photon map
  Vec3f computeRadianceWithPhotonMap(const Vec3f& wo,
                                     const SurfaceInfo& surface_info,
                                     const Primitive* hit_primitive) const {
    // get nearby photons
    const std::vector<int> photon_indices =
        photonMap.queryPhotonsInRange(surface_info.position, globalRadius);

    Vec3f Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = photonMap.getIthPhoton(photon_idx);
      const Vec3f f = hit_primitive->evaluateBxDF(
          wo, photon.wi, surface_info, TransportDirection::FROM_CAMERA);
      Lo += f * photon.throughput;
    }
    Lo /= Vec3f(nPhotons * PI * globalRadius * globalRadius);

    return Lo;
  }

  // compute in-scattering radiance with volume photon map
  Vec3f computeRadianceWithVolumePhotonMap(const Vec3f& wo, const Vec3f& pos,
                                           const Medium* medium) const {
    // get nearby photons
    const std::vector<int> photon_indices =
        volumePhotonMap.queryPhotonsInRange(pos, globalRadiusVolume);

    Vec3f Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = volumePhotonMap.getIthPhoton(photon_idx);
      const Vec3f f = medium->evalPhaseFunction(wo, photon.wi);
      Lo += f * photon.throughput;
    }
    Lo /= Vec3f(nPhotons * UNIT_SPHERE_VOLUME * globalRadiusVolume *
                globalRadiusVolume * globalRadiusVolume);

    return Lo;
  }

  // sample initial ray from light and compute initial throughput
  Ray sampleRayFromLight(const Scene& scene, Sampler& sampler) {
    // sample light
    float light_choose_pdf;
    const std::shared_ptr<Light> light =
        scene.sampleLight(sampler, light_choose_pdf);

    // sample point on light
    float light_pos_pdf;
    const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

    // sample direction on light
    float light_dir_pdf;
    const Vec3f dir =
        light->sampleDirection(light_surf, sampler, light_dir_pdf);

    // spawn ray
    Ray ray;
    ray.origin = light_surf.position;
    ray.direction = dir;
    ray.throughput = light->Le(light_surf, dir) /
                     (light_choose_pdf * light_pos_pdf * light_dir_pdf) *
                     std::abs(dot(dir, light_surf.shadingNormal));

    return ray;
  }

  // photon tracing and build photon map
  void buildPhotonMap(const Scene& scene,
                      std::vector<std::unique_ptr<Sampler>>& samplers) {
    // photon tracing
    std::vector<Photon> photons;
    std::vector<Photon> photons_volume;

#pragma omp parallel for
    for (uint32_t i = 0; i < nPhotons; ++i) {
      auto& sampler_per_thread = *samplers[omp_get_thread_num()];

      // sample initial ray from light and set initial throughput
      Ray ray = sampleRayFromLight(scene, sampler_per_thread);

      // trace photons
      uint32_t depth = 0;
      while (depth < maxDepth) {
        // invalid throughput check
        if (std::isnan(ray.throughput[0]) || std::isnan(ray.throughput[1]) ||
            std::isnan(ray.throughput[2])) {
          spdlog::error("[PPM] photon throughput is NaN");
          break;
        } else if (ray.throughput[0] < 0 || ray.throughput[1] < 0 ||
                   ray.throughput[2] < 0) {
          spdlog::error("[PPM] photon throughput is minus");
          break;
        }

        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          // russian roulette
          if (depth > 0) {
            const float russian_roulette_prob = std::min(
                std::max(ray.throughput[0],
                         std::max(ray.throughput[1], ray.throughput[2])),
                1.0f);
            if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
              break;
            }
            ray.throughput /= Vec3f(russian_roulette_prob);
          }

          // sample medium
          bool is_scattered = false;
          if (ray.hasMedium()) {
            const Medium* medium = ray.getCurrentMedium();

            Vec3f pos;
            Vec3f dir;
            Vec3f throughput_medium;
            is_scattered = medium->sampleMedium(ray, info.t, sampler_per_thread,
                                                pos, dir, throughput_medium);

            // add photon to the volume photon map when scattering occured
            if (is_scattered) {
#pragma omp critical
              {
                photons_volume.emplace_back(ray.throughput * throughput_medium,
                                            pos, -ray.direction);
              }
            }

            // advance ray
            ray.origin = pos;
            ray.direction = dir;

            // update throughput
            ray.throughput *= throughput_medium;
          }

          bool is_reflected_or_refracted = false;
          if (!is_scattered) {
            Vec3f dir = ray.direction;
            if (info.hitPrimitive->hasSurface()) {
              const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

              // add photon whenever hitting diffuse surface
              if (bxdf_type == BxDFType::DIFFUSE) {
#pragma omp critical
                {
                  photons.emplace_back(ray.throughput,
                                       info.surfaceInfo.position,
                                       -ray.direction);
                }
              }

              // sample direction by BxDF
              float pdf_dir;
              const Vec3f f = info.hitPrimitive->sampleBxDF(
                  -ray.direction, info.surfaceInfo,
                  TransportDirection::FROM_LIGHT, sampler_per_thread, dir,
                  pdf_dir);

              // update throughput
              ray.throughput *= f *
                                cosTerm(-ray.direction, dir, info.surfaceInfo,
                                        TransportDirection::FROM_LIGHT) /
                                pdf_dir;

              is_reflected_or_refracted = true;
            }

            // update ray's medium
            updateMedium(ray, dir, info);

            // update ray
            ray.origin = info.surfaceInfo.position;
            ray.direction = dir;
          }

          // update depth
          if (is_scattered || is_reflected_or_refracted) {
            depth++;
          }
        } else {
          // photon goes to the sky
          break;
        }
      }
    }

    // build photon map
    photonMap.setPhotons(photons);
    photonMap.build();

    // build volume photon map
    volumePhotonMap.setPhotons(photons_volume);
    volumePhotonMap.build();
  }

  // compute incoming radiance with photon map
  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const {
    Ray ray = ray_in;
    ray.throughput = Vec3f(1, 1, 1);

    uint32_t depth = 0;
    while (depth < maxDepth) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (depth > 0) {
          const float russian_roulette_prob =
              std::min(std::max(ray.throughput[0],
                                std::max(ray.throughput[1], ray.throughput[2])),
                       1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          ray.throughput /= Vec3f(russian_roulette_prob);
        }

        // when directly hitting light
        if (info.hitPrimitive->hasAreaLight()) {
          return ray.throughput *
                 info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

        // if hitting diffuse surface, compute reflected radiance with photon
        // map
        if (bxdf_type == BxDFType::DIFFUSE) {
          return ray.throughput *
                 computeRadianceWithPhotonMap(-ray.direction, info.surfaceInfo,
                                              info.hitPrimitive);
        }
        // if hitting specular surface, generate next ray and continue tracing
        else if (bxdf_type == BxDFType::SPECULAR) {
          // sample direction by BxDF
          Vec3f dir;
          float pdf_dir;
          Vec3f f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
              sampler, dir, pdf_dir);

          // update throughput
          ray.throughput *= f *
                            cosTerm(-ray.direction, dir, info.surfaceInfo,
                                    TransportDirection::FROM_CAMERA) /
                            pdf_dir;

          // update ray
          ray.origin = info.surfaceInfo.position;
          ray.direction = dir;
        }
      } else {
        // ray goes out the the sky
        break;
      }

      // update depth
      if (true) {
        depth++;
      }
    }

    return Vec3f(0);
  }

 public:
  PPM(const std::shared_ptr<Camera>& camera, uint32_t nIterations,
      uint32_t nPhotons, float alpha, float initialRadius,
      uint32_t maxDepth = 100)
      : Integrator(camera),
        nIterations(nIterations),
        nPhotons(nPhotons),
        alpha(alpha),
        globalRadius(initialRadius),
        maxDepth(maxDepth),
        nEmittedPhotons(0) {}

  void render(const Scene& scene, Sampler& sampler, Image& image) override {
    // init sampler for each thread
    std::vector<std::unique_ptr<Sampler>> samplers(omp_get_max_threads());
    for (int i = 0; i < samplers.size(); ++i) {
      samplers[i] = sampler.clone();
      samplers[i]->setSeed(sampler.getSeed() * (i + 1));

      // warpup sampler
      for (int j = 0; j < 10; ++j) {
        samplers[i]->getNext1D();
      }
    }

    const uint32_t width = image.getWidth();
    const uint32_t height = image.getHeight();

    spdlog::info("[PPM] rendering...");
    for (uint32_t iteration = 0; iteration < nIterations; ++iteration) {
      spdlog::info("[PPM] iteration: {}", iteration);
      spdlog::info("[PPM] radius: {}", globalRadius);

      // clear previous photon map
      photonMap.clear();
      volumePhotonMap.clear();

      // photon tracing and build photon map
      buildPhotonMap(scene, samplers);
      nEmittedPhotons += nPhotons;

      // ray tracing from camera
#pragma omp parallel for collapse(2) schedule(dynamic, 1)
      for (uint32_t i = 0; i < height; ++i) {
        for (uint32_t j = 0; j < width; ++j) {
          auto& sampler_per_thread = *samplers[omp_get_thread_num()];

          // SSAA
          const float u =
              -(2.0f * (j + sampler_per_thread.getNext1D()) - width) / height;
          const float v =
              (2.0f * (i + sampler_per_thread.getNext1D()) - height) / height;

          Ray ray;
          float pdf;
          if (camera->sampleRay(Vec2f(u, v), sampler_per_thread, ray, pdf)) {
            // compute incoming radiance with photon map
            const Vec3f radiance =
                integrate(ray, scene, sampler_per_thread) / pdf;

            // invalid radiance check
            if (std::isnan(radiance[0]) || std::isnan(radiance[1]) ||
                std::isnan(radiance[2])) {
              spdlog::error("[PPM] radiance is NaN");
              continue;
            } else if (std::isinf(radiance[0]) || std::isinf(radiance[1]) ||
                       std::isinf(radiance[2])) {
              spdlog::error("[PPM] radiance is inf");
              continue;
            } else if (radiance[0] < 0 || radiance[1] < 0 || radiance[2] < 0) {
              spdlog::error("[PPM] radiance is minus");
              continue;
            }

            // add contribution
            image.addPixel(i, j, radiance);
          } else {
            image.setPixel(i, j, Vec3f(0));
          }
        }
      }

      // update search radius
      globalRadius =
          std::sqrt((iteration + alpha) / (iteration + 1)) * globalRadius;

      // save image at each iteration
      // Image image_copied = image;
      // image_copied /= Vec3f(iteration + 1);
      // image_copied.gammaCorrection(2.2f);
      // image_copied.writePPM("iteration_" + std::to_string(iteration) +
      // ".ppm");
    }

    // take average
    image /= Vec3f(nIterations);

    spdlog::info("[PPM] done");
  }
};

#endif