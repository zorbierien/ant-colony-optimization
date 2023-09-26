//
// Created by schne on 21.09.2023.
//

#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char** readGraphFile(char *filePath, int *bufferSize) {
    FILE *file = fopen(filePath, "r");

    // Allocate Array with 10000 potential entries of nodes
    int bufferLength = 10000;
    char **lineBuffer = (char**) malloc(bufferLength * sizeof(char*));
    for (int i = 0; i < bufferLength; i++) {
        lineBuffer[i] = (char *) calloc(13, sizeof(char));
    }

    if (file) {
        printf("%s %s\n", "Try to open file ", filePath);
        printf("%s", "Datei erfolgreich geöffnet\n");

        int index = 0;
        while (!feof(file)) {
            fread(lineBuffer[index], sizeof(char), 12, file);
            if (strncmp(lineBuffer[index], "EOF", 3) == 0) {
                break;
            }
            index++;
        }
        *bufferSize = index;

        printf("%i %s\n", index, "Zeilen erfolgreich gelesen");

        // Free unnecessary memory
        for (int i = index; i < bufferLength; i++) {
            free(lineBuffer[i]);
        }
    }


    return lineBuffer;
}

graphEntry** buildAdjacenceMatrix(char **nodeArray, int arrayLength) {
    //Init AdjacenceMatrix & read Node Data in One Loop because of performance reasons
    printf("%s\n", "Build Adjacence Matrix");

    graphEntry **adjMatrix = (graphEntry**)malloc(arrayLength * sizeof(graphEntry*));
    nodeCoords nodes[arrayLength];
    for (int i = 0; i < arrayLength; i++) {
        sscanf(nodeArray[i], "%d %f %f\n", &nodes[i].id, &nodes[i].x, &nodes[i].y);
        adjMatrix[i] = (graphEntry*)malloc(arrayLength * sizeof(graphEntry));
    }

    // Calculate Euclidian Distance between all nodes and Init Pheromones
    for (int i = 0; i < arrayLength; i++) {
        for (int j = 0; j < arrayLength; j++) {
            graphEntry entry;
            float xd = nodes[i].x - nodes[j].x;
            float yd = nodes[i].y - nodes[j].y;
            entry.cost = (int)nearbyint(sqrt((double)powf(xd, 2) + (double)powf(yd, 2)));
            //TODO Check Pheromone Value -> Laut Buch "a small positive number"
            entry.pheromone = 0.00079396585946804284;
            adjMatrix[i][j] = entry;
        }
    }

    return adjMatrix;
}

int buildGraph(char *filePath, graphEntry ***adjacenceMatrix) {
    int bufferLength;
    char **lines = readGraphFile(filePath, &bufferLength);
    graphEntry **adjMatrix = buildAdjacenceMatrix(lines, bufferLength);

    //Free Memory of Read Lines
    for (int i = 0; i < bufferLength; i++) {
        free(lines[i]);
    }
    free(lines);

    *adjacenceMatrix = adjMatrix;
    return bufferLength;

// Print Adjacence Matrix
//    for (int i = 0; i < bufferLength; i++) {
//        printf("%i: ", i+1);
//        for (int j = 0; j < bufferLength; j++) {
//            printf("%f\t", adjMatrix[i][j].cost);
//        }
//        printf("\n");
//    }
}
