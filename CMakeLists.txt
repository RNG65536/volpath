cmake_minimum_required(VERSION 3.0)
project(volpath)

find_package(openmp REQUIRED)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(sky_model_srcs
  "src/hosek/ArHosekSkyModel.cpp"
  "src/hosek/ArHosekSkyModel.h"
  "src/hosek/ArHosekSkyModelData_CIEXYZ.h"
  "src/hosek/ArHosekSkyModelData_RGB.h"
  "src/hosek/ArHosekSkyModelData_Spectral.h"
)

set(srcs
    "src/camera.h"
    "src/constants.h"
    "src/fractal.h"
    "src/framebuffer.h"
    "src/framebuffer.cpp"
    "src/light.h"
    "src/main.cpp"
    "src/medium.h"
    "src/ray.h"
    "src/scene.h"
    "src/skydome.h"
    "src/texture3d.h"
    "src/texture3d.cpp"
    "src/timer.h"
    "src/timer.cpp"
    "src/utils.h"
    "src/utils.cpp"
    "src/vector.h"
    "src/vector.cpp"
    "src/volpath_cloud.h"
    "src/volpath_multichannel.h"
    "src/volpath.h"
    "src/volume.h"
    "src/volume.cpp"
)

add_library(sky_model STATIC ${sky_model_srcs})
add_executable(volpath ${srcs})

target_link_libraries(volpath sky_model)
