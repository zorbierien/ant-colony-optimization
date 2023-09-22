//
// Created by Steve on 12.09.2023.
//

#include "aco.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>


ant* initAnts(int antCount, int nodeCount) {
    if (antCount > 0) {
        ant* antArray = (ant*)malloc(antCount * sizeof(ant));
        for (int i = 0; i < antCount; i++) {
            ant newAnt;
            newAnt.length = 0;
            newAnt.hop = 0;
            newAnt.path = (int*)calloc(nodeCount, sizeof(int));
            newAnt.tabuList = (bool*)calloc(nodeCount, sizeof(bool));
            antArray[i] = newAnt;
        }
        return antArray;
    }
    else return NULL;
}

void placeAnts(ant* antArray, int antCount, int nodeCount) {
    for (int i = 0; i < antCount; i++) {
        int randomNode = ((rand()) % (nodeCount - 1) + 1);
        antArray[i].tabuList[randomNode-1] = true;
        antArray[i].currentNode = randomNode;
        antArray[i].path[0] = randomNode;
        antArray[i].hop++;
    }
}

int choosePath(ant singleAnt, graphEntry **adjMatrix, int adjMatrixLength) {
    double alpha = 1.0;
    double beta = 0.5;
    double overallPathSum = 0;
    double probabilities[adjMatrixLength];

    // Calculate sum of Pheromone * Visibilty of all possible ways
    for (int i = 0; i < adjMatrixLength; i++) {
        if (singleAnt.tabuList[i]) continue;
        else {
            overallPathSum += pow(adjMatrix[singleAnt.currentNode-1][i].pheromone, alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][i].cost, beta);
        }
    }

    // Happens if Ant is at the last node
    if (overallPathSum == 0) return singleAnt.path[0];

    // Calculate probabilities for each path
    for (int i = 0; i < adjMatrixLength; i++) {
        if (singleAnt.tabuList[i]) probabilities[i] = 0;
        else {
            probabilities[i] = (pow(adjMatrix[singleAnt.currentNode-1][i].pheromone, alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][i].cost, beta)) / overallPathSum;
        }
    }

    // Random  number between 0 and 1
    double randomNumber = (double)rand() / (double)RAND_MAX;
    double sum = 0;
    int index = 0;
    do {
        sum += probabilities[index];
        index++;
    }
    while (sum < randomNumber);
    return index;
}

void moveAnts(ant *ants, graphEntry **adjMatrix, int antCount, int adjMatrixLength) {
    for (int i = 0; i < antCount; i++) {
        int toNode = choosePath(ants[i], adjMatrix, adjMatrixLength);
        ants[i].tabuList[toNode-1] = true;
        ants[i].length += adjMatrix[ants[i].currentNode-1][toNode-1].cost;
        ants[i].currentNode = toNode;
        ants[i].path[ants[i].hop] = toNode;
        ants[i].hop++;
    }
}

int findShortestPath(ant *ants, int antCount, int **path) {
    int min = INFINITY;
    for (int i = 0; i < antCount; i++) {
        if (ants[i].length < min) min = ants[i].length;
    }
    return min;
}

void updatePheromoneLevel(graphEntry **adjacenceMatrix, int adjacenceMatrixLength, ant* ants, int antCount) {
    double ro = 0.5;
    double Q = 100;
    double deltaPheromones[adjacenceMatrixLength][adjacenceMatrixLength];
    for (int i = 0; i < adjacenceMatrixLength; i++) {
        for (int j = 0; j < adjacenceMatrixLength; j++) {
            deltaPheromones[i][j] = 0.0;
        }
    }

    for (int i = 0; i < antCount; i++) {
        double pheromoneAdd = Q / ants[i].length;
        for (int j = 0; j < adjacenceMatrixLength - 1; j++) {
            deltaPheromones[ants[i].path[j]-1][ants[i].path[j+1]-1] += pheromoneAdd;
        }
    }

    for (int i = 0; i < adjacenceMatrixLength; i++) {
        for (int j = 0; j < adjacenceMatrixLength; j++) {
            adjacenceMatrix[i][j].pheromone = (1.0-ro) * adjacenceMatrix[i][j].pheromone + deltaPheromones[i][j];
        }
    }
}

void resetAnts(ant *ants, int antCount, int nodeCount) {
    for (int i = 0; i < antCount; i++) {
        ants[i].hop = 0;
        ants[i].length = 0;
        ants[i].currentNode = 0;
        for (int j = 0; j < nodeCount; j++) {
            ants[i].tabuList[j] = false;
            ants[i].path[j] = 0;
        }
    }
}

int antColonyOptimize(char *filePath, int **path, int cycles, int numAnts) {
    time_t t;
    srand((unsigned) time(&t));
    graphEntry **adjacenceMatrix = NULL;
    int adjacenceMatrixLength = 0;
    adjacenceMatrixLength = buildGraph("../a280.tsp", &adjacenceMatrix);
    if (numAnts == 0) numAnts = adjacenceMatrixLength;
    ant* ants = initAnts(numAnts, adjacenceMatrixLength);

    int pathLength = INFINITY;
    for (int it = 0; it < cycles; it++) {
        printf("Starting cycle %d\n", it + 1);
        placeAnts(ants, numAnts, adjacenceMatrixLength);
        for (int i = 0; i < adjacenceMatrixLength; i++) {
            moveAnts(ants, adjacenceMatrix, numAnts, adjacenceMatrixLength);
        }

        int *singlePathTraverse;
        int singlePathLength = findShortestPath(ants, numAnts, &singlePathTraverse);
        if (singlePathLength < pathLength) {
            pathLength = singlePathLength;
            *path = singlePathTraverse;
        }

        updatePheromoneLevel(adjacenceMatrix, adjacenceMatrixLength, ants, numAnts);
        resetAnts(ants, numAnts, adjacenceMatrixLength);
    }

    for (int i = 0; i < adjacenceMatrixLength; i++) {
        free(adjacenceMatrix[i]);
    }
    free(adjacenceMatrix);

    for (int i = 0; i < numAnts; i++) {
        free(ants[i].path);
        free(ants[i].tabuList);
    }
    free(ants);

    return pathLength;
}
