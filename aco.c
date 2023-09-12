//
// Created by Steve on 12.09.2023.
//

#include "aco.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

colony* init_colony() {
    colony *antColony = malloc(sizeof(colony));
    antColony->cycle = 0;
    antColony->antGraph = malloc(sizeof (graph));
};

void read_graph_from_file(char* filePath, colony *antColony) {
    FILE *file = fopen(filePath, "r");
    char buffer[12];

    if (file != NULL) {
        while(!feof(file)) {
            fgets(buffer, 11, file);
            int nums[3];
            sscanf(buffer, "%d %d %d", &nums[0], &nums[1], &nums[2]);
            for (int i = 0; i < 3; i++) {
                printf("%d\n", nums[i]);
            }
        }
    }
}