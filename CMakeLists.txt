cmake_minimum_required(VERSION 3.22)
project(ACO_SEQ C)

set(CMAKE_C_STANDARD 23)

add_executable(ACO_SEQ main.c
        aco.c
        aco.h
        graph.c
        graph.h
)
target_link_libraries(ACO_SEQ m)
