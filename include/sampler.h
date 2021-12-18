#ifndef _SAMPLER_H
#define _SAMPLER_H
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "core.h"

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
typedef struct {
  uint64_t state;
  uint64_t inc;
} pcg32_random_t;

inline uint32_t pcg32_random_r(pcg32_random_t* rng) {
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// random number generator
class RNG {
 private:
  pcg32_random_t state;

 public:
  RNG() {
    state.state = 1;
    state.inc = 1;
  }
  RNG(uint64_t seed) {
    state.state = seed;
    state.inc = 1;
  }

  uint64_t getSeed() const { return state.state; }
  void setSeed(uint64_t seed) { state.state = seed; }

  float getNext() {
    constexpr float divider = 1.0f / std::numeric_limits<uint32_t>::max();
    return pcg32_random_r(&state) * divider;
  }
};

// sampler interface
class Sampler {
 protected:
  RNG rng;

 public:
  Sampler() {}

  Sampler(uint64_t seed) : rng(seed) {}

  uint64_t getSeed() const { return rng.getSeed(); }
  void setSeed(uint64_t seed) { rng.setSeed(seed); }

  virtual std::unique_ptr<Sampler> clone() const = 0;
  virtual float getNext1D() = 0;
  virtual Vec2f getNext2D() = 0;
};

// uniform distribution sampler
class UniformSampler : public Sampler {
 public:
  UniformSampler() : Sampler() {}
  UniformSampler(uint64_t seed) : Sampler(seed) {}

  std::unique_ptr<Sampler> clone() const override {
    return std::make_unique<UniformSampler>();
  }

  float getNext1D() override { return rng.getNext(); }
  Vec2f getNext2D() override { return Vec2f(rng.getNext(), rng.getNext()); }
};

// sample direction in the hemisphere
// its pdf is propotional to cosine
inline Vec3f sampleCosineHemisphere(const Vec2f& uv, float& pdf) {
  const float theta =
      0.5f * std::acos(std::clamp(1.0f - 2.0f * uv[0], -1.0f, 1.0f));
  const float phi = PI_MUL_2 * uv[1];
  const float cosTheta = std::cos(theta);
  pdf = PI_INV * cosTheta;
  return sphericalToCartesian(theta, phi);
}

// sample point on the disk
Vec2f sampleDisk(const Vec2f& uv, float R, float& pdf) {
  const float r = R * std::sqrt(std::max(uv[0], 0.0f));
  const float theta = PI_MUL_2 * uv[1];
  pdf = 1.0f / (R * R) * PI_INV;
  return Vec2f(r * std::cos(theta), r * std::sin(theta));
}

// sample value from 1D discrete empirical distribution
class DiscreteEmpiricalDistribution1D {
 private:
  std::vector<float> cdf;
  std::vector<float> pdf;

 public:
  DiscreteEmpiricalDistribution1D(const float* values, unsigned int N) {
    // sum f
    float sum = 0;
    for (std::size_t i = 0; i < N; ++i) {
      sum += values[i];
    }

    // compute cdf
    cdf.resize(N + 1);
    cdf[0] = 0;
    for (std::size_t i = 1; i < N + 1; ++i) {
      cdf[i] = cdf[i - 1] + values[i - 1] / sum;
    }

    // compute pdf
    pdf.resize(N);
    for (std::size_t i = 0; i < N; ++i) {
      pdf[i] = cdf[i + 1] - cdf[i];
    }
  }

  DiscreteEmpiricalDistribution1D(const std::vector<float>& values)
      : DiscreteEmpiricalDistribution1D(values.data(), values.size()) {}

  uint32_t sample(float u, float& pdf) const {
    // inverse cdf
    int x = std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    if (x == 0) {
      x++;
    }

    // compute pdf
    pdf = cdf[x] - cdf[x - 1];

    // NOTE: cdf's index is +1 from values
    return x - 1;
  }

  float getPDF(uint32_t i) const { return pdf[i]; }
};

#endif