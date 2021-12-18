#ifndef _RAY_H
#define _RAY_H
#include <stack>

#include "core.h"

// forward declaration
class Medium;

class Ray {
 public:
  Vec3f origin;
  Vec3f direction;
  Vec3f throughput;
  std::stack<const Medium*> mediums;

  static constexpr float tmin = RAY_EPS;
  float tmax = std::numeric_limits<float>::max();

  Ray() {}
  Ray(const Vec3f& origin, const Vec3f& direction)
      : origin(origin), direction(direction) {}

  Vec3f operator()(float t) const { return origin + t * direction; }

  bool hasMedium() const { return !mediums.empty(); }

  const Medium* getCurrentMedium() const { return mediums.top(); }

  void pushMedium(const Medium* medium) { mediums.push(medium); }

  void popMedium() {
    if (hasMedium()) {
      mediums.pop();
    }
  }
};

#endif