#ifndef _IMAGE_H
#define _IMAGE_H
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "core.h"

class Image {
 private:
  uint32_t width;
  uint32_t height;
  std::vector<float> pixels;

 private:
  uint32_t getIndex(uint32_t i, uint32_t j) const {
    return 3 * j + 3 * width * i;
  }

 public:
  Image(uint32_t width, uint32_t height) : width(width), height(height) {
    pixels.resize(3 * width * height);
  }

  uint32_t getWidth() const { return width; }
  uint32_t getHeight() const { return height; }

  Vec3f getPixel(uint32_t i, uint32_t j) const {
    const uint32_t idx = getIndex(i, j);
    return Vec3f(pixels[idx], pixels[idx + 1], pixels[idx + 2]);
  }

  void addPixel(uint32_t i, uint32_t j, const Vec3f& rgb) {
    const uint32_t idx = getIndex(i, j);
    pixels[idx] += rgb[0];
    pixels[idx + 1] += rgb[1];
    pixels[idx + 2] += rgb[2];
  }

  void setPixel(uint32_t i, uint32_t j, const Vec3f& rgb) {
    const uint32_t idx = getIndex(i, j);
    pixels[idx] = rgb[0];
    pixels[idx + 1] = rgb[1];
    pixels[idx + 2] = rgb[2];
  }

  Image& operator*=(const Vec3f& rgb) {
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        const Vec3f c = getPixel(i, j);
        setPixel(i, j, c * rgb);
      }
    }
    return *this;
  }

  Image& operator/=(const Vec3f& rgb) {
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        const Vec3f c = getPixel(i, j);
        setPixel(i, j, c / rgb);
      }
    }
    return *this;
  }

  void gammaCorrection(const float gamma) {
    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        Vec3f c = getPixel(i, j);

        c[0] = std::pow(c[0], 1.0f / gamma);
        c[1] = std::pow(c[1], 1.0f / gamma);
        c[2] = std::pow(c[2], 1.0f / gamma);

        setPixel(i, j, c);
      }
    }
  }

  void writePPM(const std::string& filename) {
    std::ofstream file(filename);

    file << "P3" << std::endl;
    file << width << " " << height << std::endl;
    file << "255" << std::endl;

    for (int i = 0; i < height; ++i) {
      for (int j = 0; j < width; ++j) {
        const Vec3f rgb = getPixel(i, j);
        const uint32_t R =
            std::clamp(static_cast<uint32_t>(255.0f * rgb[0]), 0u, 255u);
        const uint32_t G =
            std::clamp(static_cast<uint32_t>(255.0f * rgb[1]), 0u, 255u);
        const uint32_t B =
            std::clamp(static_cast<uint32_t>(255.0f * rgb[2]), 0u, 255u);
        file << R << " " << G << " " << B << std::endl;
      }
    }

    file.close();
  }
};

#endif