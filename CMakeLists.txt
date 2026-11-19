cmake_minimum_required(VERSION 3.20)
project(LAMBDA VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(LAMBDA_FETCH_DEPS   "Download Eigen/yaml-cpp if not found on the system" ON)
option(LAMBDA_BUILD_TESTS  "Build tests" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(Dependencies)

add_library(lambda_core
  src/io/Config.cpp
  src/io/Output.cpp
  src/io/FieldWriter.cpp
  src/post/Eigenfunction.cpp
  src/solver/RadialBasis.cpp
  src/solver/Solver.cpp)
target_include_directories(lambda_core PUBLIC src)
target_link_libraries(lambda_core PUBLIC lambda::deps)
target_compile_definitions(lambda_core PUBLIC LAMBDA_VERSION="${PROJECT_VERSION}")
if(MSVC)
  target_compile_options(lambda_core PUBLIC /W3 /utf-8 /bigobj)
  target_compile_definitions(lambda_core PUBLIC _USE_MATH_DEFINES NOMINMAX)
else()
  target_compile_options(lambda_core PRIVATE -Wall -Wextra)
endif()

add_executable(lambdaDrop src/main.cpp)
target_link_libraries(lambdaDrop PRIVATE lambda_core)
install(TARGETS lambdaDrop RUNTIME DESTINATION bin)

if(LAMBDA_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
