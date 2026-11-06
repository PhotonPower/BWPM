cmake_minimum_required(VERSION 3.24)
project(BWPM VERSION 0.1.0 DESCRIPTION "Fast bidirectional vector wave propagation method" LANGUAGES NONE)

option(BWPM_BUILD_CORE "Build the compiled module bwpm._core" ON)
option(BWPM_ENABLE_CUDA "Build the CUDA backend (RTX 6000 Ada, sm_89), stage B4d" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT BWPM_BUILD_CORE)
  message(STATUS "BWPM: core disabled, pure Python reference only (no compiler needed)")
  return()
endif()
enable_language(CXX)

add_library(bwpm_core STATIC cpp/src/interface.cpp)
target_include_directories(bwpm_core PUBLIC cpp/include)
set_target_properties(bwpm_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(MSVC)
  target_compile_options(bwpm_core PRIVATE /W4)
else()
  target_compile_options(bwpm_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(SKBUILD)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core cpp/python/bindings.cpp)
  target_link_libraries(_core PRIVATE bwpm_core)
  install(TARGETS _core LIBRARY DESTINATION bwpm)
endif()

message(STATUS "BWPM ${PROJECT_VERSION}: core=${BWPM_BUILD_CORE} CUDA=${BWPM_ENABLE_CUDA}")
