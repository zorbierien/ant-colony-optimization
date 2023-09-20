//
// Created by Steve on 12.09.2023.
//

#ifndef ACO_SEQ_ACO_H
#define ACO_SEQ_ACO_H

typedef struct {
    float cost;
    double pheromone;
} graphEntry;

typedef struct {
    int id;
    float x;
    float y;
} nodeCoords;

void buildGraph(char * filePath);
char** readGraphFile(char *filePath, int *bufferSize);
graphEntry** buildAdjacenceMatrix(char **nodeArray, int arrayLength);

#endif //ACO_SEQ_ACO_H
