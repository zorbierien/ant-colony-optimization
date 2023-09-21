//
// Created by Steve on 12.09.2023.
//

#include "aco.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


ant* initAnts(int antCount) {
    if (antCount > 0) {
        ant* antArray = (ant*)malloc(antCount * sizeof(ant));
        for (int i = 0; i < antCount; i++) {
            ant newAnt;
            newAnt.length = 0;
            antArray[i] = newAnt;
        }
        return antArray;
    }
    else return NULL;
}

void placeAnts(ant* antArray, int antCount, int nodeCount) {
    time_t t;
    srand((unsigned) time(&t));
    for (int i = 0; i < antCount; i++) {
        antArray[i].tabuList = (int*) malloc(nodeCount * sizeof(int));
        antArray[i].tabuList[0] = (rand() + 1) % 280;
    }
}

void moveAnts(ant* ants) {

}

int findShortestPath(ant* ants, int* path) {

}

void updatePheromoneLevel(graphEntry **adjacenceMatrix) {

}

int antColonyOptimize(char *filePath, int *path, int cycles, int numAnts) {
    graphEntry **adjacenceMatrix = NULL;
    int adjacenceMatrixLength = 0;
    adjacenceMatrixLength = buildGraph("../a280.tsp", adjacenceMatrix);
    if (numAnts == 0) numAnts = adjacenceMatrixLength;
    ant* ants = initAnts(numAnts);
    placeAnts(ants, numAnts, adjacenceMatrixLength);
}
