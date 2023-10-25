#include <stdio.h>
#include "aco.h"
#include <sys/time.h>

int main() {
    int *path = NULL;
    int length = 0;

    // Time measurement start
    struct timeval beforeTime;
    gettimeofday(&beforeTime, NULL);
    unsigned long long beforeTimestamp = (unsigned long long)beforeTime.tv_sec * 1000000 + beforeTime.tv_usec;

    length = antColonyOptimize("../k100.tsp", &path, 1000, 100);

    // Time measurement end
    struct timeval afterTime;
    gettimeofday(&afterTime, NULL);
    unsigned long long afterTimestamp = (unsigned long long)afterTime.tv_sec * 1000000 + afterTime.tv_usec;

    printf("Shortest Path: %i\n", length);
    printf("Duration: %llu\n", afterTimestamp - beforeTimestamp);
    return 0;
}
