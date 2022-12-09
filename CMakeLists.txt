cmake_minimum_required(VERSION 3.21)
project(dgemm)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(dgemm
        dgemm.cpp
        )
