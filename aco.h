//
// Created by Steve on 12.09.2023.
//

#ifndef ACO_SEQ_ACO_H
#define ACO_SEQ_ACO_H
#include "graph.h"
#include <stdbool.h>

typedef struct {
    bool *tabuList;
    int *path;
    int hop;
    int currentNode;
    int length;
} ant;

ant* initAnts(int antCount, int nodeCount);
void placeAnts(ant* antArray, int antCount, int nodeCount);
void moveAnts(ant *ants, graphEntry **adjMatrix, int antCount, int adjMatrixLength);
int findShortestPath(ant *ants, int antCount, int **path);
void updatePheromoneLevel(graphEntry **adjacenceMatrix, const int *path, int pathArrayLength, long pathLength);
int antColonyOptimize(char *filePath, int **path, int cycles, int numAnts);

#endif //ACO_SEQ_ACO_H
