#include <iostream>
#include <chrono>
#include <vector>

void two_loops(int imax, int jmax) {
    for (int i = 0; i < imax; i++) {
        for (int j = 0; j < jmax; j++) {
            // Your operation
//			if(i == 0 && j == 0) std::cout << "-----------------------------------" << j << endl;
            std::cout << "i: " << i << " ; " << "j: " << j << endl;
        }
    }
}

void one_loop(int imax, int jmax) { // Works well for different objects
    for (int k = 0; k < imax * jmax; k++) {
        int i = k / jmax;
        int j = k % jmax;
        //if(j <= i) continue; // This makes j = i+1;
        // Your operation
//		if(k == 0) std::cout << "-----------------------------------" << j << endl;
        std::cout << "i: " << i << " ; " << "j: " << j << endl;

    }
}

int timingtest() {

    std::vector<int> vec(3, 0);
    std::vector<int> vec2(5, 0);

    auto start = std::chrono::high_resolution_clock::now();
    two_loops(vec.size(),vec2.size());
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Two loops took: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
              << " ns\n";

    start = std::chrono::high_resolution_clock::now();
    one_loop(vec.size(),vec2.size());
    end = std::chrono::high_resolution_clock::now();
    std::cout << "One loop took: "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()
              << " ns\n";

    return 0;
}
