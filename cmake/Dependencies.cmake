# Dependencies of LAMBDA.
#   Eigen 3.4 (header only) and yaml-cpp >= 0.7. No BLAS/LAPACK is needed.
# System packages are used first; with LAMBDA_FETCH_DEPS=ON missing Eigen/yaml-cpp are downloaded
# as release archives (checked by SHA-256; git is not needed).
include(FetchContent)

add_library(lambda_deps INTERFACE)
add_library(lambda::deps ALIAS lambda_deps)

# ---- Eigen ----
find_package(Eigen3 3.4 CONFIG QUIET)
if(NOT Eigen3_FOUND)
  if(NOT LAMBDA_FETCH_DEPS)
    message(FATAL_ERROR "Eigen 3.4 not found (set Eigen3_DIR or LAMBDA_FETCH_DEPS=ON)")
  endif()
  message(STATUS "Eigen: fetching 3.4.0")
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_SUBDIR _none_)          # header only: do not run Eigen's own CMake
  FetchContent_MakeAvailable(eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${eigen_SOURCE_DIR})
endif()
target_link_libraries(lambda_deps INTERFACE Eigen3::Eigen)

# ---- yaml-cpp ----
find_package(yaml-cpp 0.7 CONFIG QUIET)
if(NOT yaml-cpp_FOUND)
  if(NOT LAMBDA_FETCH_DEPS)
    message(FATAL_ERROR "yaml-cpp not found (set yaml-cpp_DIR or LAMBDA_FETCH_DEPS=ON)")
  endif()
  message(STATUS "yaml-cpp: fetching 0.8.0")
  set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
  set(YAML_CPP_INSTALL OFF CACHE BOOL "" FORCE)
  set(YAML_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(yaml-cpp
    URL https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz
    URL_HASH SHA256=fbe74bbdcee21d656715688706da3c8becfd946d92cd44705cc6098bb23b3a16
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  # yaml-cpp 0.8.0 declares cmake_minimum_required(VERSION 3.4), which CMake >= 4.0 rejects
  set(CMAKE_POLICY_VERSION_MINIMUM 3.10)
  FetchContent_MakeAvailable(yaml-cpp)
  unset(CMAKE_POLICY_VERSION_MINIMUM)
endif()
if(TARGET yaml-cpp::yaml-cpp)
  target_link_libraries(lambda_deps INTERFACE yaml-cpp::yaml-cpp)
else()
  target_link_libraries(lambda_deps INTERFACE yaml-cpp)   # older config files
endif()

