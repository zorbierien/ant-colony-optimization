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
#include <immintrin.h>

struct mmasParams {
    float ro;
    float alpha;
    float beta;
    int minimumTourLength;
};

struct mmasParams params = {.ro = 0.5f, .alpha = 1.0f, .beta= 0.5f, .minimumTourLength = 21282};

ant* initAnts(int antCount, int nodeCount) {
    if (antCount > 0) {
        ant* antArray = (ant*)malloc(antCount * sizeof(ant));
        for (int i = 0; i < antCount; i++) {
            ant newAnt;
            newAnt.length = 0;
            newAnt.hop = 0;
            newAnt.path = (int*)calloc(nodeCount+1, sizeof(int));
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
    double overallPathSum = 0;
    double probabilities[adjMatrixLength];

    // Calculate sum of Pheromone * Visibilty of all possible ways
    // OPTIMIZED: Einzelne Heuristik hier schon ion probabilities Array eingefügt, sodass anschließend nur noch eine Division notwendig
    __m256d pathSum = _mm256_set1_pd(0);
    __m256d probabilityVectors[adjMatrixLength / 4 + 1];
    double vecArray[4];
    for (int i = 0; i < adjMatrixLength / 4; i++) {
        if (singleAnt.tabuList[4*i]) vecArray[0] = 0;
        else {
            vecArray[0] = pow(adjMatrix[singleAnt.currentNode-1][4*i].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][4*i].cost, params.beta);
        }
        if (singleAnt.tabuList[4*i+1]) vecArray[1] = 0;
        else {
            vecArray[1] = pow(adjMatrix[singleAnt.currentNode-1][4*i+1].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][4*i+1].cost, params.beta);
        }
        if (singleAnt.tabuList[4*i+2]) vecArray[2] = 0;
        else {
            vecArray[2] = pow(adjMatrix[singleAnt.currentNode-1][4*i+2].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][4*i+2].cost, params.beta);
        }
        if (singleAnt.tabuList[4*i+3]) vecArray[3] = 0;
        else {
            vecArray[3] = pow(adjMatrix[singleAnt.currentNode-1][4*i+3].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][4*i+3].cost, params.beta);
        }
        __m256d vector = _mm256_set_pd(vecArray[3], vecArray[2], vecArray[1], vecArray[0]);
        probabilityVectors[i] = vector;
        pathSum = _mm256_add_pd(pathSum, vector);
    }

    if (adjMatrixLength % 4 != 0) {
        int todo = adjMatrixLength % 4;
        for (int i = 0; i < todo; i++) {
            vecArray[i] = pow(adjMatrix[singleAnt.currentNode-1][i].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][i].cost, params.beta);
        }
        for (int i = todo; i < 4; i++) {
            vecArray[i] = 0;
        }
        __m256d vector = _mm256_set_pd(vecArray[3], vecArray[2], vecArray[1], vecArray[0]);
        pathSum = _mm256_add_pd(pathSum, vector);
        probabilityVectors[adjMatrixLength / 4] = vector;
    }
    else {
        probabilityVectors[adjMatrixLength / 4] = _mm256_set1_pd(0);
    }

    double pathSumValues[4];
    _mm256_storeu_pd((double*)&pathSumValues, pathSum);
    overallPathSum = pathSumValues[0] + pathSumValues[1] + pathSumValues[2] + pathSumValues[3];

    // Happens if Ant is at the last node
    if (overallPathSum == 0) return singleAnt.path[0];

    // Calculate probabilities for each path
    //OPTIMIZED: probability wird hier nur noch durch overallPathSum geteilt, statt vollkommen neu berechnet
    pathSum = _mm256_set1_pd(overallPathSum);
    for (int i = 0; i <= adjMatrixLength / 4; i++) {
        _mm256_storeu_pd((double*)&probabilities[4*i], _mm256_div_pd(probabilityVectors[i], pathSum));
    }

    // Random  number between 0 and 1
    double randomNumber = (double)rand() / (double)RAND_MAX;
    double sum = 0;
    int index = -1;
    do {
        sum += probabilities[index+1];
        index++;
    }
    while (sum < randomNumber);

    return index+1;
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
        if (ants[i].length < min) {
            min = ants[i].length;
            *path = ants[i].path;
        }
    }
    return min;
}

// For Symmetric TSP-Problem
void updatePheromoneLevel(graphEntry **adjacenceMatrix, const int adjacenceMatrixLength, const int *path, int pathArrayLength, long pathLength) {
    // Sollte bekannt sein
    double pheromoneMin = 1.0f / (params.ro * (float)params.minimumTourLength * (float)pathArrayLength * (float)pathArrayLength);
    double pheromoneMax = 1.0f / (params.ro * (float)params.minimumTourLength);

    //OPTIMIZED: Innere Schleife nur von i+1 bis j um Symmetrie der Matrix zu nutzen
    for (int i = 0; i < adjacenceMatrixLength; ++i) {
        adjacenceMatrix[i][i].pheromone = (1.0-params.ro) * adjacenceMatrix[i][i].pheromone;
        for (int j = i + 1; j < adjacenceMatrixLength; ++j) {
            double newPheromone = (1.0-params.ro) * adjacenceMatrix[i][j].pheromone;
            adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = newPheromone;
            if (adjacenceMatrix[i][j].pheromone < pheromoneMin) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMin;
            if (adjacenceMatrix[i][j].pheromone > pheromoneMax) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMax;
        }
    }

    double pheromoneAdd = 1 / (double)pathLength;
    for (int i = 0; i < pathArrayLength - 1; ++i) {
        adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone = adjacenceMatrix[path[i+1]-1][path[i]-1].pheromone = adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone + pheromoneAdd;
        if (adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone > pheromoneMax) adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone = adjacenceMatrix[path[i+1]-1][path[i]-1].pheromone = pheromoneMax;
        if (adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone < pheromoneMin) adjacenceMatrix[path[i]-1][path[i+1]-1].pheromone = adjacenceMatrix[path[i+1]-1][path[i]-1].pheromone = pheromoneMin;
    }
}

void resetAnts(ant *ants, int antCount, int nodeCount) {
    for (int i = 0; i < antCount; i++) {
        ants[i].hop = 0;
        ants[i].length = 0;
        ants[i].currentNode = 0;
        for (int j = 0; j < nodeCount; j++) {
            ants[i].tabuList[j] = false;
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

        updatePheromoneLevel(adjacenceMatrix, adjacenceMatrixLength,singlePathTraverse, adjacenceMatrixLength+1, singlePathLength);
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
