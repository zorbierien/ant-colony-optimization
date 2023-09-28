cmake_minimum_required(VERSION 3.22)
project(ACO_SEQ C)

set(CMAKE_C_STANDARD 23)
set(OPTIMIZE_FLAGS "-O3 -Ofast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native ${OPTIMIZE_FLAGS}")

add_executable(ACO_SEQ main.c
        aco.c
        aco.h
        graph.c
        graph.h
)
target_link_libraries(ACO_SEQ ${CMAKE_SOURCE_DIR}/lib/amd-libm/lib/libamdlibm.a)
target_link_libraries(ACO_SEQ ${CMAKE_SOURCE_DIR}/lib/amd-libm/lib/libamdlibmfast.a)
target_link_libraries(ACO_SEQ m)
