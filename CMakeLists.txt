cmake_minimum_required(VERSION 3.24)
project(dgemm)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(dgemm
        dgemm.cpp
        MMIOReader.cpp MMIOReader.hpp Matrix.hpp Vector.hpp mmio.c mmio.h)
