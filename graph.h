//
// Created by schne on 21.09.2023.
//

#ifndef ACO_SEQ_GRAPH_H
#define ACO_SEQ_GRAPH_H

typedef struct {
    int cost;
    double pheromone;
} graphEntry;

typedef struct {
    int id;
    float x;
    float y;
} nodeCoords;

int buildGraph(char * filePath, graphEntry ***adjacenceMatrix);
char** readGraphFile(char *filePath, int *bufferSize);
graphEntry** buildAdjacenceMatrix(char **nodeArray, int arrayLength);

#endif //ACO_SEQ_GRAPH_H
