#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H
#include <cmath>
#include <memory>

#include "light.h"
#include "material.h"
#include "medium.h"
#include "triangle.h"

// primitive provides an abstraction layer of the objects in the scene
class Primitive {
 private:
  const Triangle* triangle;
  const BxDF* bxdf;
  const Medium* medium;
  const Light* areaLight;

 public:
  Primitive(const Triangle* triangle, const BxDF* bxdf,
            const Medium* medium = nullptr, const Light* areaLight = nullptr)
      : triangle(triangle), bxdf(bxdf), medium(medium), areaLight(areaLight) {}

  bool hasSurface() const { return bxdf != nullptr; }
  bool hasMedium() const { return medium != nullptr; }
  bool hasAreaLight() const { return areaLight != nullptr; }

  // return emission
  Vec3f Le(const SurfaceInfo& surfInfo, const Vec3f& dir) const {
    return areaLight->Le(surfInfo, dir);
  }

  const Medium* getMedium() const { return medium; }

  BxDFType getBxDFType() const { return bxdf->getType(); }

  Vec3f evaluateBxDF(const Vec3f& wo, const Vec3f& wi,
                     const SurfaceInfo& surfInfo,
                     const TransportDirection& mode) const {
    // world to local transform
    const Vec3f wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);
    const Vec3f wi_l =
        worldToLocal(wi, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);

    return bxdf->evaluate(wo_l, wi_l, mode);
  }

  // sample direction by BxDF
  // its pdf is propotional to the shape of BxDF
  Vec3f sampleBxDF(const Vec3f& wo, const SurfaceInfo& surfInfo,
                   const TransportDirection& mode, Sampler& sampler, Vec3f& wi,
                   float& pdf) const {
    // world to local transform
    const Vec3f wo_l =
        worldToLocal(wo, surfInfo.dpdu, surfInfo.shadingNormal, surfInfo.dpdv);

    // sample direction in tangent space
    Vec3f wi_l;
    const Vec3f f = bxdf->sampleDirection(wo_l, mode, sampler, wi_l, pdf);

    // local to world transform
    wi = localToWorld(wi_l, surfInfo.dpdu, surfInfo.shadingNormal,
                      surfInfo.dpdv);

    return f;
  }
};

#endif