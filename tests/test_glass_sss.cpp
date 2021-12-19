#include "camera.h"
#include "integrator.h"
#include "medium.h"
#include "scene.h"

int main() {
  const uint32_t width = 512;
  const uint32_t height = 512;
  const uint32_t n_iterations = 100;
  const uint32_t n_photons = 500000;
  const float alpha = 3.0f / 4.0f;
  const float initial_radius_surface = 0.1f;
  const float initial_radius_volume = 5.0f * initial_radius_surface;
  const uint32_t max_depth = 100;

  Image image(width, height);

  const Vec3f camera_pos = Vec3f(1, 3, 6);
  const Vec3f camera_lookat = Vec3f(-0.254774, 0.653768, -0.172298);
  const Vec3f camera_forward = normalize(camera_lookat - camera_pos);
  const float FOV = 0.25 * PI;

  const auto camera =
      std::make_shared<PinholeCamera>(camera_pos, camera_forward, FOV);

  // build scene
  Scene scene;
  scene.loadObj("bunny-test.obj");
  scene.build();

  // render
  UniformSampler sampler;
  PPM integrator(camera, n_iterations, n_photons, alpha, initial_radius_surface,
                 initial_radius_volume, max_depth);
  // PathTracing integrator(camera, 10000, 101);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}