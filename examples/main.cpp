#include "camera.h"
#include "integrator.h"
#include "scene.h"

int main() {
  const uint32_t width = 512;
  const uint32_t height = 512;
  const uint32_t n_iterations = 1000;
  const uint32_t n_photons = 200000;
  const float alpha = 3.0f / 4.0f;
  const float initial_radius_surface = 0.01f;
  const float initial_radius_volume = 5.0f * initial_radius_surface;
  const uint32_t max_depth = 100;

  Image image(width, height);

  const Vec3f camera_position = Vec3f(0, 1, 6);
  const Vec3f lookat = Vec3f(0, 1, 0);
  const Vec3f forward = normalize(lookat - camera_position);
  const auto camera =
      std::make_shared<PinholeCamera>(camera_position, forward, 0.25 * PI);
  // const auto camera = std::make_shared<ThinLensCamera>(camera_position,
  // forward,
  //                                                      0.38f * PI, 64.0f);
  // camera->focus(Vec3f(0, 1, 0.5));

  // build scene
  Scene scene;
  scene.loadObj(
      "models/CornellBox-Water-Small-Light-Covered-Cube-Floor-No-Backwall.obj");
  scene.build();

  // render
  UniformSampler sampler;
  // PathTracing integrator(camera, 100);
  // integrator.render(scene, sampler, image);
  PPM integrator(camera, n_iterations, n_photons, alpha, initial_radius_surface,
                 initial_radius_volume, max_depth);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}