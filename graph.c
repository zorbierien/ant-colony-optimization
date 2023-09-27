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
    int bufferLength = 8;
    int strlen = 20;
    char **lineBuffer = (char**) malloc(bufferLength * sizeof(char*));
    for (int i = 0; i < bufferLength; i++) {
        lineBuffer[i] = (char *) calloc(strlen, sizeof(char));
    }

    if (file) {
        printf("%s %s\n", "Try to open file ", filePath);
        printf("%s", "Datei erfolgreich geöffnet\n");

        int index = 0;
        while (!feof(file)) {
            // OPTIMIZED Allocate more buffer if needed
            if (bufferLength < index + 1) {
                lineBuffer = (char**)realloc(lineBuffer, bufferLength * 2 * sizeof (char*));
                if (lineBuffer == NULL) exit(1);
                for (int i = bufferLength; i < bufferLength * 2; i++) {
                    lineBuffer[i] = (char *) calloc(strlen, sizeof(char));
                }
                bufferLength *= 2;
            }
            getline(&lineBuffer[index], (size_t*)&strlen, file);
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
    // OPTIMIZED: Use the Symmetry of the matrix
    for (int i = 0; i < arrayLength; i++) {
        graphEntry doubleEntry;
        doubleEntry.cost = 0;
        doubleEntry.pheromone = 100;
        adjMatrix[i][i] = doubleEntry;
        for (int j = i + 1; j < arrayLength; j++) {
            graphEntry entry;
            float xd = nodes[i].x - nodes[j].x;
            float yd = nodes[i].y - nodes[j].y;
            entry.cost = (int)nearbyint(sqrt((double)powf(xd, 2) + (double)powf(yd, 2)));
            entry.pheromone = 100;
            adjMatrix[i][j] = entry;

            graphEntry symmetricEntry;
            symmetricEntry.cost = entry.cost;
            symmetricEntry.pheromone = entry.pheromone;
            adjMatrix[j][i] = symmetricEntry;
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
}
