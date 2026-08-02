cmake_minimum_required(VERSION 3.24)
project(Lab_1)

set(CMAKE_C_STANDARD 99)

add_executable(Lab_1 main.c include/matrix.h src/matrix.c src/algorithms.c include/matrix_slice.h src/matrix_slice.c)
