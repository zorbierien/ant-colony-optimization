#include "aco.h"
#include "graph.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

struct mmasParams {
    float ro;
    float alpha;
    float beta;
    // Should be known or a good approximation
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
    prob probabilities[adjMatrixLength];

    // Calculate sum of Pheromone * Visibilty of all possible ways
    double singleHeuristic;
    for (int i = 0; i < adjMatrixLength; i++) {
        if (singleAnt.tabuList[i]) continue;
        else {
            singleHeuristic = pow(adjMatrix[singleAnt.currentNode-1][i].pheromone, params.alpha) * pow(1.0 / adjMatrix[singleAnt.currentNode-1][i].cost, params.beta);
            overallPathSum += singleHeuristic;
            probabilities[i].probability = singleHeuristic;
        }
    }

    // Happens if Ant is at the last node
    if (overallPathSum == 0) return singleAnt.path[0];

    // Calculate probabilities for each path
    for (int i = 0; i < adjMatrixLength; i++) {
        if (singleAnt.tabuList[i]) {
            probabilities[i].probability = 0;
            probabilities[i].node = i + 1;
        }
        else {
            probabilities[i].probability /= overallPathSum;
            probabilities[i].node = i + 1;
        }
    }

    // Random  number between 0 and 1
    double randomNumber = (double)rand() / (double)RAND_MAX;
    double sum = 0;
    int index = -1;
    // Sum up probabilities to random number to choose path
    do {
        sum += probabilities[index+1].probability;
        index++;
    }
    while (sum < randomNumber);

    return probabilities[index].node;
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
