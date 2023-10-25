#include "aco.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <immintrin.h>
#include <limits.h>
#define AMD_LIBM_VEC_EXPERIMENTAL
#include "lib/amd-libm/include/amdlibm.h"
#include "lib/amd-libm/include/amdlibm_vec.h"

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
    __m256d pathSum = _mm256_set1_pd(0);
    __m256d probabilityVectors[adjMatrixLength / 4 + 1];
    __m256d alphaVector = _mm256_set1_pd(params.alpha);
    __m256d betaVector = _mm256_set1_pd(params.beta);
    for (int i = 0; i < adjMatrixLength / 4; i++) {
        double pheromoneArray[4];
        double costArray[4];
        __m256d pheromoneVector;
        __m256d costVector;
        if (singleAnt.tabuList[4*i]) {
            costArray[0] = INT_MAX;
            pheromoneArray[0] = 0;
        }
        else {
            costArray[0] = adjMatrix[singleAnt.currentNode-1][4*i].cost;
            pheromoneArray[0] = adjMatrix[singleAnt.currentNode-1][4*i].pheromone;
        }
        if (singleAnt.tabuList[4*i+1]) {
            costArray[1] = INT_MAX;
            pheromoneArray[1] = 0;
        }
        else {
            costArray[1] = adjMatrix[singleAnt.currentNode-1][4*i+1].cost;
            pheromoneArray[1] = adjMatrix[singleAnt.currentNode-1][4*i+1].pheromone;
        }
        if (singleAnt.tabuList[4*i+2]) {
            costArray[2] = INT_MAX;
            pheromoneArray[2] = 0;
        }
        else {
            costArray[2] = adjMatrix[singleAnt.currentNode-1][4*i+2].cost;
            pheromoneArray[2] = adjMatrix[singleAnt.currentNode-1][4*i+2].pheromone;
        }
        if (singleAnt.tabuList[4*i+3]) {
            costArray[3] = INT_MAX;
            pheromoneArray[3] = 0;
        }
        else {
            costArray[3] = adjMatrix[singleAnt.currentNode-1][4*i+3].cost;
            pheromoneArray[3] = adjMatrix[singleAnt.currentNode-1][4*i+3].pheromone;
        }
        costVector = _mm256_loadu_pd((double*)&costArray);
        __m256d visibiltyVector = _mm256_div_pd(_mm256_set1_pd(1), costVector);
        pheromoneVector = _mm256_loadu_pd((double*)&pheromoneArray);
        probabilityVectors[i] = _mm256_mul_pd(amd_vrd4_pow(pheromoneVector, alphaVector), amd_vrd4_pow(visibiltyVector, betaVector));
        pathSum = _mm256_add_pd(pathSum, probabilityVectors[i]);
    }

    if (adjMatrixLength % 4 != 0) {
        int todo = adjMatrixLength % 4;
        for (int i = 0; i < todo; i++) {
            probabilities[adjMatrixLength - 4 + i] = amd_pow(adjMatrix[singleAnt.currentNode-1][i].pheromone, params.alpha) * amd_pow(1.0 / adjMatrix[singleAnt.currentNode-1][i].cost, params.beta);
            overallPathSum += probabilities[adjMatrixLength - 4 + i];
        }
        for (int i = todo; i < 4; i++) {
            probabilities[adjMatrixLength - 4 + i] = 0;
        }
    }

    double pathSumValues[4];
    _mm256_storeu_pd((double*)&pathSumValues, pathSum);
    overallPathSum = pathSumValues[0] + pathSumValues[1] + pathSumValues[2] + pathSumValues[3];

    // Happens if Ant is at the last node
    if (overallPathSum == 0) return singleAnt.path[0];

    // Calculate probabilities for each path
    pathSum = _mm256_set1_pd(overallPathSum);
    for (int i = 0; i < adjMatrixLength / 4; i++) {
        _mm256_storeu_pd((double*)&probabilities[4*i], _mm256_div_pd(probabilityVectors[i], pathSum));
    }

    for (int i = 0; i < adjMatrixLength % 4; i++) {
        probabilities[adjMatrixLength / 4 + i] /= overallPathSum;
    }

    // Random  number between 0 and 1
    double randomNumber = (double)rand() / (double)RAND_MAX;
    double sum = 0;
    int index = -1;
    // Sum up probabilities to random number to choose path
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
    int min = INT_MAX;
    for (int i = 0; i < antCount; i++) {
        if (ants[i].length < min) {
            min = ants[i].length;
            *path = ants[i].path;
        }
    }
    return min;
}

void updatePheromoneLevel(graphEntry **adjacenceMatrix, const int adjacenceMatrixLength, const int *path, int pathArrayLength, long pathLength) {
    double pheromoneMin = 1.0f / (params.ro * (float)params.minimumTourLength * (float)pathArrayLength * (float)pathArrayLength);
    double pheromoneMax = 1.0f / (params.ro * (float)params.minimumTourLength);

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
    adjacenceMatrixLength = buildGraph(filePath, &adjacenceMatrix);
    if (numAnts == 0) numAnts = adjacenceMatrixLength;
    ant* ants = initAnts(numAnts, adjacenceMatrixLength);

    int pathLength = INT_MAX;
    for (int it = 0; it < cycles; it++) {
        placeAnts(ants, numAnts, adjacenceMatrixLength);
        for (int i = 0; i < adjacenceMatrixLength; i++) {
            moveAnts(ants, adjacenceMatrix, numAnts, adjacenceMatrixLength);
        }

        int *singlePathTraverse;
        int singlePathLength = findShortestPath(ants, numAnts, &singlePathTraverse);
        if (singlePathLength < pathLength) {
            pathLength = singlePathLength;
            *path = singlePathTraverse;
            printf("New best route with length %d found in iteration %d\n", singlePathLength, it);
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
