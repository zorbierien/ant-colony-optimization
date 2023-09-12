//
// Created by Steve on 12.09.2023.
//

#ifndef ACO_SEQ_ACO_H
#define ACO_SEQ_ACO_H

typedef struct {
    int *tabuList;
} ant;

typedef struct {
  int numNodes;
  int *adjList
} graph;

typedef struct {
    graph *antGraph;
    int cycle;
} colony;

colony* init_colony();
void read_graph_from_file(char* filePath, colony *antColony);

#endif //ACO_SEQ_ACO_H
