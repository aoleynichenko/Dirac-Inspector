#
# DIRAC Inspector
# 2024 A. Oleynichenko
#

cmake_minimum_required(VERSION 3.19)
enable_language(C)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(C_STANDARD_REQUIRED ON)
project(dirac_inspector)

add_executable(dirac_inspector.x
        src/main.c
        src/mdprop.c
        src/mrconee.c
        src/libunf.c
        src/mdcint.c
        src/mdcint.h
)

target_link_libraries(dirac_inspector.x -lm)



