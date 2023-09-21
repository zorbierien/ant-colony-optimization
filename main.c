#include <stdio.h>
#include "aco.h"

int main() {
    int *path = NULL;
    int length = 0;
    length = antColonyOptimize("../a280.tsp", path, 1, 280);
    return 0;
}
