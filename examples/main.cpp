#include "camera.h"
#include "integrator.h"
#include "scene.h"

int main() {
  const uint32_t width = 512;
  const uint32_t height = 512;
  const uint32_t n_iterations = 10;
  const uint32_t n_photons = 100000;
  const float alpha = 3.0f / 4.0f;
  const float initial_radius = 0.01f;
  const uint32_t max_depth = 100;

  Image image(width, height);

  const auto camera = std::make_shared<ThinLensCamera>(
      Vec3f(0, 1, 6), Vec3f(0, 0, -1), 0.38f * PI, 32.0f);
  camera->focus(Vec3f(-0.2496, -0.001, 0.6));

  // build scene
  Scene scene;
  scene.loadObj("cornellbox-water2.obj");
  scene.build();

  // render
  UniformSampler sampler;
  // PathTracing integrator(camera, 1000);
  // integrator.render(scene, sampler, image);
  PPM integrator(camera, n_iterations, n_photons, alpha, initial_radius,
                 max_depth);
  integrator.render(scene, sampler, image);

  // gamma correction
  image.gammaCorrection(2.2f);

  // output image
  image.writePPM("output.ppm");

  return 0;
}