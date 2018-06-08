//
// Created by Michael Heuer on 08.06.18.
//
#include <stdio.h>
#include <omp.h>
#include <vector>
#include <iostream>

int main(int argc, char *argv[]) {
    int omp_get_thread_num();
# pragma omp parallel
    {
        printf("Thread rank: %d\n", omp_get_thread_num());
    }

    int n = 100;
    std::vector<int> numbers;
    numbers.reserve(n);

    double wtime = omp_get_wtime();
# pragma omp parallel for shared(numbers)
    for (int i = 0; i < n; ++i) {
        numbers[i] = i;
    }
    wtime = omp_get_wtime()-wtime;
    std::cout << "\"elapsed wall-clock time: " << wtime << " seconds" << std::endl;

    for (int i = 0; i < n; ++i) {
        std::cout << numbers[i] << std::endl;
    }

}