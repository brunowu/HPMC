cmake_minimum_required (VERSION 3.8)

project (SparseMatrix LANGUAGES C)

set(CMAKE_C_FLAGS_RELEASE "-O3")

list (APPEND CMAKE_MODULE_PATH "${SparseMatrix_SOURCE_DIR}/cmake-modules")

set(BLA_VENDOR Intel10_64lp)
find_package(BLAS)
find_package(OpenMP)

find_package(CUDAToolkit)

if(MKL_FOUND)
	message("-- MKL has been detected")
endif()

if(OpenMP_FOUND)
	message("-- System has support for OpenMP")
else()
	message("-- System has no support for OpenMP")
endif()

if(CUDAToolkit_FOUND)
	message("-- CUDA has been detected")
endif()

add_subdirectory(SparseMatrix)

add_subdirectory(hands-on)

add_subdirectory(homework-1)

option(BUILD_HOMEWORK2 "Build the homework2" OFF)

if(BUILD_HOMEWORK2)
	message("======>> Build with homework 2")
   if(BLAS_FOUND)
   add_subdirectory(homework-2)
   endif()
endif()


if(CUDAToolkit_FOUND)
  add_subdirectory(bonus)
endif()
