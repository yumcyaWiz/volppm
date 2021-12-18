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
              (2.0f * (j + sampler_per_pixel->getNext1D()) - width) / height;
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

// implementation of path tracing
// NOTE: for reference purpose
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
    Vec3f throughput(1, 1, 1);

    for (uint32_t k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // russian roulette
        if (k > 0) {
          const float russian_roulette_prob = std::min(
              std::max(throughput[0], std::max(throughput[1], throughput[2])),
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob) {
            break;
          }
          throughput /= russian_roulette_prob;
        }

        // Le
        if (info.hitPrimitive->hasAreaLight()) {
          radiance += throughput *
                      info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        // sample direction by BxDF
        Vec3f dir;
        float pdf_dir;
        Vec3f f = info.hitPrimitive->sampleBxDF(
            -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
            sampler, dir, pdf_dir);

        // update throughput and ray
        throughput *= f *
                      cosTerm(-ray.direction, dir, info.surfaceInfo,
                              TransportDirection::FROM_CAMERA) /
                      pdf_dir;
        ray = Ray(info.surfaceInfo.position, dir);
      } else {
        break;
      }
    }

    return radiance;
  }
};

// this implementation is based on modified version of original SPPM
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
  // global search radius for radiance estimation
  float globalRadius;

  PhotonMap photonMap;

  // compute reflected radiance with photon map
  Vec3f computeRadianceWithPhotonMap(const Vec3f& wo,
                                     const IntersectInfo& info) const {
    // get nearby photons
    const std::vector<int> photon_indices =
        photonMap.queryPhotonsInRange(info.surfaceInfo.position, globalRadius);

    Vec3f Lo;
    for (const int photon_idx : photon_indices) {
      const Photon& photon = photonMap.getIthPhoton(photon_idx);
      const Vec3f f = info.hitPrimitive->evaluateBxDF(
          wo, photon.wi, info.surfaceInfo, TransportDirection::FROM_CAMERA);
      Lo += f * photon.throughput;
    }
    Lo /= (nPhotons * PI * globalRadius * globalRadius);

    return Lo;
  }

  // sample initial ray from light and compute initial throughput
  Ray sampleRayFromLight(const Scene& scene, Sampler& sampler,
                         Vec3f& throughput) {
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
    Ray ray(light_surf.position, dir);
    throughput = light->Le(light_surf, dir) /
                 (light_choose_pdf * light_pos_pdf * light_dir_pdf) *
                 std::abs(dot(dir, light_surf.shadingNormal));

    return ray;
  }

  // photon tracing and build photon map
  void buildPhotonMap(const Scene& scene,
                      std::vector<std::unique_ptr<Sampler>>& samplers) {
    // photon tracing
    std::vector<Photon> photons;

    // spdlog::info("[PPMAPA] tracing photons...");
#pragma omp parallel for
    for (uint32_t i = 0; i < nPhotons; ++i) {
      auto& sampler_per_thread = *samplers[omp_get_thread_num()];

      // sample initial ray from light and set initial throughput
      Vec3f throughput;
      Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

      // trace photons
      // whener hitting diffuse surface, add photon to the photon array
      // recursively tracing photon with russian roulette
      for (uint32_t k = 0; k < maxDepth; ++k) {
        if (std::isnan(throughput[0]) || std::isnan(throughput[1]) ||
            std::isnan(throughput[2])) {
          spdlog::error("[PPM] photon throughput is NaN");
          break;
        } else if (throughput[0] < 0 || throughput[1] < 0 ||
                   throughput[2] < 0) {
          spdlog::error("[PPM] photon throughput is minus");
          break;
        }

        IntersectInfo info;
        if (scene.intersect(ray, info)) {
          const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
          if (bxdf_type == BxDFType::DIFFUSE) {
            // TODO: remove lock to get more speed
#pragma omp critical
            {
              photons.emplace_back(throughput, info.surfaceInfo.position,
                                   -ray.direction);
            }
          }

          // russian roulette
          if (k > 0) {
            const float russian_roulette_prob = std::min(
                std::max(throughput[0], std::max(throughput[1], throughput[2])),
                1.0f);
            if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
              break;
            }
            throughput /= russian_roulette_prob;
          }

          // sample direction by BxDF
          Vec3f dir;
          float pdf_dir;
          const Vec3f f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, TransportDirection::FROM_LIGHT,
              sampler_per_thread, dir, pdf_dir);

          // update throughput and ray
          throughput *= f *
                        cosTerm(-ray.direction, dir, info.surfaceInfo,
                                TransportDirection::FROM_LIGHT) /
                        pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        } else {
          // photon goes to the sky
          break;
        }
      }
    }
    // spdlog::info("[PPMAPA] done");

    // build photon map
    // spdlog::info("[PPMAPA] building photon map...");
    photonMap.setPhotons(photons);
    photonMap.build();
    // spdlog::info("[PPMAPA] done");
  }

  // compute incoming radiance with photon map
  Vec3f integrate(const Ray& ray_in, const Scene& scene,
                  Sampler& sampler) const {
    Ray ray = ray_in;
    Vec3f throughput(1, 1, 1);

    for (uint32_t k = 0; k < maxDepth; ++k) {
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // when directly hitting light
        if (info.hitPrimitive->hasAreaLight()) {
          return throughput *
                 info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

        // if hitting diffuse surface, compute reflected radiance with photon
        // map
        if (bxdf_type == BxDFType::DIFFUSE) {
          return throughput *
                 computeRadianceWithPhotonMap(-ray.direction, info);
        }
        // if hitting specular surface, generate next ray and continue tracing
        else if (bxdf_type == BxDFType::SPECULAR) {
          // sample direction by BxDF
          Vec3f dir;
          float pdf_dir;
          Vec3f f = info.hitPrimitive->sampleBxDF(
              -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
              sampler, dir, pdf_dir);

          // update throughput and ray
          throughput *= f *
                        cosTerm(-ray.direction, dir, info.surfaceInfo,
                                TransportDirection::FROM_CAMERA) /
                        pdf_dir;
          ray = Ray(info.surfaceInfo.position, dir);
        }
      } else {
        // ray goes out the the sky
        break;
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

      // photon tracing and build photon map
      // spdlog::info("[PPMAPA] photon tracing pass...");
      buildPhotonMap(scene, samplers);
      nEmittedPhotons += nPhotons;
      // spdlog::info("[PPMAPA] done");

      // eye tracing
      // spdlog::info("[PPMAPA] eye tracing pass...");
#pragma omp parallel for collapse(2) schedule(dynamic, 1)
      for (uint32_t i = 0; i < height; ++i) {
        for (uint32_t j = 0; j < width; ++j) {
          auto& sampler_per_thread = *samplers[omp_get_thread_num()];

          // SSAA
          const float u =
              (2.0f * (j + sampler_per_thread.getNext1D()) - width) / height;
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
      // spdlog::info("[SPPM] done");

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