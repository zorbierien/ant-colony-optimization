//
// Created by Steve on 12.09.2023.
//

#ifndef ACO_SEQ_ACO_H
#define ACO_SEQ_ACO_H
#include "graph.h"

typedef struct {
    int *tabuList;
    int length;
} ant;

ant* initAnts(int antCount);
void placeAnts(ant* antArray, int antCount, int nodeCount);
void moveAnts(ant* ants);
int findShortestPath(ant* ants, int* path);
void updatePheromoneLevel(graphEntry **adjacenceMatrix);
int antColonyOptimize(char *filePath, int *path, int cycles, int numAnts);

#endif //ACO_SEQ_ACO_H
