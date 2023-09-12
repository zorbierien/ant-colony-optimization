#include <stdio.h>
#include "aco.h"

int main() {
    colony *antColony = init_colony();
    read_graph_from_file("a280.tsp", antColony);
    return 0;
}
