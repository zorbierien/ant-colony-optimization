#include <stdio.h>
#include "aco.h"
#include <sys/time.h>

int main() {
    int *path = NULL;
    int length = 0;

    struct timeval beforeTime;
    gettimeofday(&beforeTime, NULL);
    unsigned long long beforeTimestamp = (unsigned long long)beforeTime.tv_sec * 1000000 + beforeTime.tv_usec;

    length = antColonyOptimize("../a280.tsp", &path, 10, 280);

    struct timeval afterTime;
    gettimeofday(&beforeTime, NULL);
    unsigned long long afterTimestamp = (unsigned long long)beforeTime.tv_sec * 1000000 + beforeTime.tv_usec;

    printf("Shortest Path: %i\n", length);
    printf("Duration: %llu\n", afterTimestamp - beforeTimestamp);
    return 0;
}
