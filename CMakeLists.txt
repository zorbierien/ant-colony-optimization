cmake_minimum_required(VERSION 3.26)
project(ACO_SEQ C)

set(CMAKE_C_STANDARD 23)

add_executable(ACO_SEQ main.c
        aco.c
        aco.h)
