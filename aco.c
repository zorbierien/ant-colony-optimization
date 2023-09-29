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
#include <limits.h>
#include <string.h>

struct mmasParams {
    float ro;
    float alpha;
    float beta;
    int minimumTourLength;
};

struct mmasParams params = {.ro = 0.5f, .alpha = 1.0f, .beta= 0.5f, .minimumTourLength = 22000};

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
    __m256d probabilityVectors[adjMatrixLength / 4];
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
        __m256d vector = _mm256_loadu_pd((double*)&vecArray);
        probabilityVectors[i] = vector;
        pathSum = _mm256_add_pd(pathSum, vector);
    }

    if (adjMatrixLength % 4 != 0) {
        int todo = adjMatrixLength % 4;
        for (int i = 0; i < todo; i++) {
            int arrayIndex = (adjMatrixLength / 4) * 4 + i;
            if (singleAnt.tabuList[arrayIndex]) {
                probabilities[arrayIndex] = 0;
                continue;
            }
            probabilities[arrayIndex] = pow(adjMatrix[singleAnt.currentNode-1][arrayIndex].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][arrayIndex].cost, params.beta);
            overallPathSum += probabilities[arrayIndex];
        }
    }

    double pathSumValues[4];
    _mm256_storeu_pd((double*)&pathSumValues, pathSum);
    overallPathSum += pathSumValues[0] + pathSumValues[1] + pathSumValues[2] + pathSumValues[3];

    // Happens if Ant is at the last node
    if (overallPathSum == 0) return singleAnt.path[0];

    // Calculate probabilities for each path
    //OPTIMIZED: probability wird hier nur noch durch overallPathSum geteilt, statt vollkommen neu berechnet
    pathSum = _mm256_set1_pd(overallPathSum);
    for (int i = 0; i < adjMatrixLength / 4; i++) {
        __m256d normProbability = _mm256_div_pd(probabilityVectors[i], pathSum);
        _mm256_storeu_pd((double*)&probabilities[4*i], normProbability);
    }

    for (int i = 0; i < adjMatrixLength % 4; i++) {
        probabilities[(adjMatrixLength / 4) * 4 + i] /= overallPathSum;
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
    int min = INT_MAX;
    for (int i = 0; i < antCount; i++) {
        if (ants[i].length < min) {
            min = ants[i].length;
            *path = ants[i].path;
        }
    }
    return min;
}

// For Symmetric TSP-Problem
//TODO noch checken
void updatePheromoneLevel(graphEntry **adjacenceMatrix, const int adjacenceMatrixLength, const int *path, int pathArrayLength, long pathLength) {
    double pheromoneMin = 1.0f / (params.ro * (float)params.minimumTourLength * (float)pathArrayLength * (float)pathArrayLength);
    double pheromoneMax = 1.0f / (params.ro * (float)params.minimumTourLength);

    __m256d evaporationVector = _mm256_set1_pd(1.0-params.ro);
    for (int i = 0; i < adjacenceMatrixLength; ++i) {
        if (i % 4 == 0 && i <= adjacenceMatrixLength - 4) {
            __m256d adjacencePheromoneVector = _mm256_set_pd(adjacenceMatrix[i+3][i+3].pheromone, adjacenceMatrix[i+2][i+2].pheromone, adjacenceMatrix[i+1][i+1].pheromone, adjacenceMatrix[i][i].pheromone);
            adjacencePheromoneVector = _mm256_mul_pd(evaporationVector, adjacencePheromoneVector);
            __m256d pheromoneMinCmp = _mm256_cmp_pd(adjacencePheromoneVector, _mm256_set1_pd(pheromoneMin), _CMP_LT_OQ);
            __m256d pheromoneMaxCmp = _mm256_cmp_pd(adjacencePheromoneVector, _mm256_set1_pd(pheromoneMax), _CMP_GT_OQ);
            double pheromoneMinArray[4];
            _mm256_storeu_pd((double*)&pheromoneMinArray, pheromoneMinCmp);
            double pheromoneMaxArray[4];
            _mm256_storeu_pd((double*)&pheromoneMaxArray, pheromoneMaxCmp);
            double newPheromoneArray[4];
            _mm256_storeu_pd((double*)&newPheromoneArray, adjacencePheromoneVector);
            if (pheromoneMinArray[0] != 0) adjacenceMatrix[i+3][i+3].pheromone = pheromoneMin;
            else adjacenceMatrix[i+3][i+3].pheromone = newPheromoneArray[0];
            if (pheromoneMinArray[1] != 0) adjacenceMatrix[i+2][i+2].pheromone = pheromoneMin;
            else adjacenceMatrix[i+2][i+2].pheromone = newPheromoneArray[1];
            if (pheromoneMinArray[2] != 0) adjacenceMatrix[i+1][i+1].pheromone = pheromoneMin;
            else adjacenceMatrix[i+1][i+1].pheromone = newPheromoneArray[2];
            if (pheromoneMinArray[3] != 0) adjacenceMatrix[i][i].pheromone = pheromoneMin;
            else adjacenceMatrix[i][i].pheromone = newPheromoneArray[3];
            if (pheromoneMaxArray[0] != 0) adjacenceMatrix[i+3][i+3].pheromone = pheromoneMax;
            if (pheromoneMaxArray[1] != 0) adjacenceMatrix[i+2][i+2].pheromone = pheromoneMax;
            if (pheromoneMaxArray[2] != 0) adjacenceMatrix[i+1][i+1].pheromone = pheromoneMax;
            if (pheromoneMaxArray[3] != 0) adjacenceMatrix[i][i].pheromone = pheromoneMax;
        }
        else if (i >= adjacenceMatrixLength - 4 + adjacenceMatrixLength % 4 && adjacenceMatrixLength % 4 != 0) {
            adjacenceMatrix[i][i].pheromone = adjacenceMatrix[i][i].pheromone * (1.0-params.ro);
            if (adjacenceMatrix[i][i].pheromone < pheromoneMin) adjacenceMatrix[i][i].pheromone = pheromoneMin;
            if (adjacenceMatrix[i][i].pheromone > pheromoneMax) adjacenceMatrix[i][i].pheromone = pheromoneMax;
        }

        for (int j = i + 1; j < adjacenceMatrixLength; ++j) {
            if((j - (i+1)) % 4 == 0 && j <= adjacenceMatrixLength - 4) {
                __m256d adjacencePheromoneVector = _mm256_set_pd(adjacenceMatrix[i][j+3].pheromone, adjacenceMatrix[i][j+2].pheromone, adjacenceMatrix[i][j+1].pheromone, adjacenceMatrix[i][j].pheromone);
                adjacencePheromoneVector = _mm256_mul_pd(evaporationVector, adjacencePheromoneVector);
                __m256d pheromoneMinCmp = _mm256_cmp_pd(adjacencePheromoneVector, _mm256_set1_pd(pheromoneMin), _CMP_LT_OQ);
                __m256d pheromoneMaxCmp = _mm256_cmp_pd(adjacencePheromoneVector, _mm256_set1_pd(pheromoneMax), _CMP_GT_OQ);
                double pheromoneMinArray[4];
                _mm256_storeu_pd((double*)&pheromoneMinArray, pheromoneMinCmp);
                double pheromoneMaxArray[4];
                _mm256_storeu_pd((double*)&pheromoneMaxArray, pheromoneMaxCmp);
                double newPheromoneArray[4];
                _mm256_storeu_pd((double*)&newPheromoneArray, adjacencePheromoneVector);
                if (pheromoneMinArray[0] != 0) adjacenceMatrix[i][j+3].pheromone = adjacenceMatrix[j+3][i].pheromone = pheromoneMin;
                else adjacenceMatrix[i][j+3].pheromone = adjacenceMatrix[j+3][i].pheromone = newPheromoneArray[0];
                if (pheromoneMinArray[1] != 0) adjacenceMatrix[i][j+2].pheromone = adjacenceMatrix[j+2][i].pheromone = pheromoneMin;
                else adjacenceMatrix[i][j+2].pheromone = adjacenceMatrix[j+2][i].pheromone = newPheromoneArray[1];
                if (pheromoneMinArray[2] != 0) adjacenceMatrix[i][j+1].pheromone = adjacenceMatrix[j+1][i].pheromone = pheromoneMin;
                else adjacenceMatrix[i][j+1].pheromone = adjacenceMatrix[j+1][i].pheromone =  newPheromoneArray[2];
                if (pheromoneMinArray[3] != 0) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMin;
                else adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = newPheromoneArray[3];
                if (pheromoneMaxArray[0] != 0) adjacenceMatrix[i][j+3].pheromone = adjacenceMatrix[j+3][i].pheromone = pheromoneMax;
                if (pheromoneMaxArray[1] != 0) adjacenceMatrix[i][j+2].pheromone = adjacenceMatrix[j+2][i].pheromone = pheromoneMax;
                if (pheromoneMaxArray[2] != 0) adjacenceMatrix[i][j+1].pheromone = adjacenceMatrix[j+1][i].pheromone = pheromoneMax;
                if (pheromoneMaxArray[3] != 0) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMax;
            }
            else if (j >= adjacenceMatrixLength - ((adjacenceMatrixLength - (i+1)) % 4) && (adjacenceMatrixLength - (i+1)) % 4 != 0) {
                adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = adjacenceMatrix[i][j].pheromone * (1.0-params.ro);
                if (adjacenceMatrix[i][j].pheromone < pheromoneMin) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMin;
                if (adjacenceMatrix[i][j].pheromone > pheromoneMax) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMax;
            }
        }
    }

//    for (int i = 0; i < adjacenceMatrixLength && adjacenceMatrixLength % 4 != 0; i++) {
//        for (int j = 0; j < adjacenceMatrixLength % 4; j++) {
//            adjacenceMatrix[i][adjacenceMatrixLength - 1 - j].pheromone = adjacenceMatrix[adjacenceMatrixLength - 1 - j][i].pheromone = adjacenceMatrix[i][adjacenceMatrixLength - 1 - j].pheromone * (1.0-params.ro);
//        }
//    }

    //OPTIMIZED: Innere Schleife nur von i+1 bis j um Symmetrie der Matrix zu nutzen
//    for (int i = 0; i < adjacenceMatrixLength; ++i) {
//        adjacenceMatrix[i][i].pheromone = (1.0-params.ro) * adjacenceMatrix[i][i].pheromone;
//        for (int j = i + 1; j < adjacenceMatrixLength; ++j) {
//            double newPheromone = (1.0-params.ro) * adjacenceMatrix[i][j].pheromone;
//            adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = newPheromone;
//            if (adjacenceMatrix[i][j].pheromone < pheromoneMin) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMin;
//            if (adjacenceMatrix[i][j].pheromone > pheromoneMax) adjacenceMatrix[i][j].pheromone = adjacenceMatrix[j][i].pheromone = pheromoneMax;
//        }
//    }

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

    FILE *file = fopen("../test.txt", "a");
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
            printf("Better solution %d found at iteration %d\n", pathLength, it);
        }

        updatePheromoneLevel(adjacenceMatrix, adjacenceMatrixLength,singlePathTraverse, adjacenceMatrixLength+1, singlePathLength);
        for (int i = 0; i < adjacenceMatrixLength; i++) {
            char text[200];
            sprintf(text, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", adjacenceMatrix[i][0].pheromone, adjacenceMatrix[i][1].pheromone, adjacenceMatrix[i][2].pheromone, adjacenceMatrix[i][3].pheromone, adjacenceMatrix[i][4].pheromone, adjacenceMatrix[i][5].pheromone, adjacenceMatrix[i][6].pheromone, adjacenceMatrix[i][7].pheromone, adjacenceMatrix[i][8].pheromone, adjacenceMatrix[i][9].pheromone);
            int n = (int)strlen(text);
            fwrite(&text, sizeof(char), strlen(text), file);
        }
        char newLine[1] = {"\n"};
        fwrite(&newLine, sizeof(char), 1, file);
        resetAnts(ants, numAnts, adjacenceMatrixLength);
    }
    fclose(file);
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
