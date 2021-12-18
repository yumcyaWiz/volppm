#ifndef _CAMERA_H
#define _CAMERA_H
#include <spdlog/spdlog.h>

#include <cmath>

#include "core.h"
#include "ray.h"
#include "sampler.h"

class Camera {
 protected:
  Vec3f position;
  Vec3f forward;
  Vec3f right;
  Vec3f up;

 public:
  Camera(const Vec3f& position, const Vec3f& forward)
      : position(position), forward(forward) {
    right = normalize(cross(forward, Vec3f(0, 1, 0)));
    up = normalize(cross(right, forward));

    spdlog::info("[Camera] position: ({}, {}, {})", position[0], position[1],
                 position[2]);
    spdlog::info("[Camera] forward: ({}, {}, {})", forward[0], forward[1],
                 forward[2]);
    spdlog::info("[Camera] right: ({}, {}, {})", right[0], right[1], right[2]);
    spdlog::info("[Camera] up: ({}, {}, {})", up[0], up[1], up[2]);
  }

  // sample ray emitting from the given sensor coordinate
  // NOTE: uv: [-1, -1] x [1, 1], sensor coordinate
  virtual bool sampleRay(const Vec2f& uv, Sampler& sampler, Ray& ray,
                         float& pdf) const = 0;
};

// pinhole camera
class PinholeCamera : public Camera {
 private:
  float FOV;
  float focalLength;

 public:
  PinholeCamera(const Vec3f& position, const Vec3f& forward,
                float FOV = 0.5f * PI)
      : Camera(position, forward) {
    // compute focal length from FOV
    focalLength = 1.0f / std::tan(0.5f * FOV);
    spdlog::info("[PinholeCamera] focalLength: {}", focalLength);
  }

  bool sampleRay(const Vec2f& uv, Sampler& sampler, Ray& ray,
                 float& pdf) const override {
    const Vec3f pinholePos = position + focalLength * forward;
    const Vec3f sensorPos = position + uv[0] * right + uv[1] * up;
    ray = Ray(sensorPos, normalize(pinholePos - sensorPos));
    pdf = 1.0f;
    return true;
  }
};

// thin lens camera
class ThinLensCamera : public Camera {
 private:
  float focalLength;
  float lensRadius;
  // distance from the sensor to the lens
  float a;
  // distance from the lens to the object plane
  float b;

 public:
  ThinLensCamera(const Vec3f& position, const Vec3f& forward,
                 float FOV = 0.5f * PI, float fNumber = 8.0f)
      : Camera(position, forward) {
    // compute focal length from FOV
    focalLength = 1.0f / std::tan(0.5f * FOV);
    spdlog::info("[ThinLensCamera] focalLength: {}", focalLength);

    // compute lens radius from F-number
    lensRadius = 2.0f * focalLength / fNumber;
    spdlog::info("[ThinLensCamera] lensRadius: {}", lensRadius);

    // init a, b
    // focus at inf
    b = 10000.0f;
    a = 1.0f / (1.0f / focalLength - 1.0f / b);
    spdlog::info("[ThinLensCamera] a: {}", a);
    spdlog::info("[ThinLensCamera] b: {}", b);
  }

  // focus at specified position
  void focus(const Vec3f& p) {
    b = dot(p - position, forward) - a;
    a = 1.0f / (1.0f / focalLength - 1.0f / b);

    spdlog::info("[ThinLensCamera] focusing at ({}, {}, {})", p[0], p[1], p[2]);
    spdlog::info("[ThinLensCamera] a: {}", a);
    spdlog::info("[ThinLensCamera] b: {}", b);
  }

  bool sampleRay(const Vec2f& uv, Sampler& sampler, Ray& ray,
                 float& pdf) const override {
    const Vec3f sensorPos = position + uv[0] * right + uv[1] * up;
    const Vec3f lensCenter = position + a * forward;

    // sample point on lens
    float pdf_area;
    const Vec2f pLens2D = sampleDisk(sampler.getNext2D(), lensRadius, pdf_area);
    const Vec3f pLens = lensCenter + pLens2D[0] * right + pLens2D[1] * up;
    Vec3f sensorToLens = normalize(pLens - sensorPos);

    // find intersection point with object plane
    const Vec3f sensorToLensCenter = normalize(lensCenter - sensorPos);
    const Vec3f pObject =
        sensorPos +
        ((a + b) / dot(sensorToLensCenter, forward)) * sensorToLensCenter;

    ray = Ray(pLens, normalize(pObject - pLens));
    // pdf = length2(pLens - sensorPos) / dot(sensorToLens, forward) * pdf_area;
    pdf = 1.0f;
    return true;
  }
};

#endif