cmake_minimum_required(VERSION 3.20)
project(volppm LANGUAGES CXX)

option(BUILD_TESTS "build tests" OFF)

# spdlog
if(NOT TARGET spdlog)
    find_package(spdlog REQUIRED)
endif()

# OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(WARNING "failed to find OpenMP")
endif()

# Embree3
find_package(embree 3.0 REQUIRED)
if (${embree_FOUND})
  message(STATUS "Found Embree")
else()
  message(FATAL_ERROR "Could not find Embree")
endif()

# externals
add_subdirectory("externals")

# volppm
add_library(volppm INTERFACE)
target_include_directories(volppm INTERFACE "include")
target_compile_features(volppm INTERFACE cxx_std_20)
set_target_properties(volppm PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(volppm INTERFACE spdlog::spdlog)
target_link_libraries(volppm INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(volppm INTERFACE tinyobjloader)
target_link_libraries(volppm INTERFACE embree)

# examples
add_subdirectory("examples")

# tests
if(BUILD_TESTS)
  add_subdirectory("tests")
endif()