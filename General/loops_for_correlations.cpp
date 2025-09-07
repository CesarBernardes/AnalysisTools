///*** To compile use:
///*** g++ -std=c++17 loops_for_correlations.cpp -o loops_for_correlations
///*** To run use:
///*** ./loops_for_correlations 

#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <cmath>

void two_loops(int imax, int jmax) {
    for (int i = 0; i < imax; i++) {
        for (int j = i+1; j < jmax; j++) {
//            std::cout << "i: " << i << " ; " << "j: " << j << std::endl;
        }
    }
}

void one_loop(int imax, int jmax) { // Works well for different objects
    for (int k = 0; k < imax * jmax; k++) {
        int i = k / jmax;
        int j = k % jmax;
        if(j <= i) continue; // This makes j = i+1;
//        std::cout << "i: " << i << " ; " << "j: " << j << std::endl;

    }
}

void two_loops_parallel(int imax, int jmax, int thread_count) {
    std::vector<std::thread> threads;
    
    // Função que cada thread executará
    auto thread_task = [imax, jmax](int start_i, int end_i) {
        for (int i = start_i; i < end_i; i++) {
            for (int j = i+1; j < jmax; j++) {
//                std::cout << "i: " << i << " ; " << "j: " << j << std::endl;
            }
        }
    };
    
    // Distribui o trabalho entre as threads
    int chunk_size = imax / thread_count;
    for (int t = 0; t < thread_count; t++) {
        int start_i = t * chunk_size;
        int end_i = (t == thread_count - 1) ? imax : (t + 1) * chunk_size;
        threads.emplace_back(thread_task, start_i, end_i);
    }
    
    // Aguarda todas as threads terminarem
    for (auto& thread : threads) {
        thread.join();
    }
}



int main() {

    std::vector<int> vec(10000, 0);
    std::vector<int> vec2(10000, 0);
    int thread_count = std::thread::hardware_concurrency(); // Número de cores disponíveis

    // Medição do two_loops original
    auto start = std::chrono::high_resolution_clock::now();
    two_loops(vec.size(),vec2.size());
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Two loops took: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*std::pow(10,-9)
              << " s\n";

    // Medição do one_loop
    start = std::chrono::high_resolution_clock::now();
    one_loop(vec.size(),vec2.size());
    end = std::chrono::high_resolution_clock::now();
    std::cout << "One loop took: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*std::pow(10,-9)
              << " s\n";

    // Medição do two_loops paralelizado
    start = std::chrono::high_resolution_clock::now();
    two_loops_parallel(vec.size(), vec2.size(), thread_count);
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Two loops (parallel, " << thread_count << " threads) took: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*std::pow(10,-9)
              << " s\n";

    return 0;
}
