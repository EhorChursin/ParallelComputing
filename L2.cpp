
#include <iostream>
#include <fstream>
#include <vector>
#include <stdbool.h>
#include <cstdlib>
#include "time.h"
#include <omp.h>


void saveMatrixToFile(const std::vector<std::vector<int>>& matrix, const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    std::streampos fileSize = 0;

    int numRows = matrix.size();
    int numCols = (numRows > 0) ? matrix[0].size() : 0;

    for (const auto& row : matrix) {
        file.write((char*)row.data(), sizeof(int) * row.size());
        fileSize += sizeof(int) * row.size();
    }
    file.close();
    std::cout << filename << " saved, total memory: " << fileSize << " bytes, rows: " << numRows << ", columns: " << numCols << std::endl;
}

std::vector<std::vector<int>> loadMatrixFromFile(const std::string& filename, int rows, int cols) {
    std::vector<std::vector<int>> matrix(rows, std::vector<int>(cols));
    std::ifstream file(filename, std::ios::binary);
    for (auto& row : matrix) {
        file.read((char*)row.data(), sizeof(int) * cols);
    }
    file.close();
    return matrix;
}

void printMatrix(const std::vector<std::vector<int>>& matrix) {
    std::cout << std::endl;
    for (const auto& row : matrix) {
        for (int num : row) {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<int>> matrixMultOpenMP(const std::vector<std::vector<int>>& m1, const std::vector<std::vector<int>>& m2) {

    int rows1 = m1.size();
    int cols1 = m1[0].size();
    int cols2 = m2[0].size();

    std::vector<std::vector<int>> result(rows1, std::vector<int>(cols2, 0));

    int threads;
    omp_set_num_threads(4);
    #pragma omp parallel shared(threads)
        {
            threads = omp_get_num_threads();

    #pragma for mp parallel shared(threads)

        for (int i = 0; i < rows1; i++) {
            for (int j = 0; j < cols2; j++) {
                for (int k = 0; k < cols1; k++) {
                    result[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
    }

    return result;
}

void saveResult(const std::vector<std::vector<int>>& matrix, const double& time, const int& r1, const int& c1, const int& c2) {

    std::string filename = "res-omp.txt";
    std::ofstream file(filename);

    double volume = r1 * c1 * c2;

    for (const auto& row : matrix) {
        for (int num : row) {
            file << num << " ";
        }
        file << "\n";
    }

    file << "\n";
    file << "Calculation time on 4 threads: " << time << " seconds." << std::endl;
    file << "Number of multiplication operations: " << volume << std::endl;
    file << "Number of saves to memory: " << volume << std::endl;
    file << "Number of downloads from memory: " << 2 * volume << std::endl;
    file.close();
}

int main() {

    srand(time(nullptr));

    bool SaveToFile = false;

    std::vector<int> s1 = { 800, 200 };
    std::vector<int> s2 = { 200, 400 };

    std::string m1_filename = "m1.bin";
    std::string m2_filename = "m2.bin";

    if (SaveToFile) {

        std::vector<std::vector<int>> m1(s1[0], std::vector<int>(s1[1]));
        std::vector<std::vector<int>> m2(s2[0], std::vector<int>(s2[1]));

        for (int i = 0; i < s1[0]; i++) {
            for (int j = 0; j < s1[1]; j++) {
                m1[i][j] = rand() % 100;
            }
        }

        for (int i = 0; i < s2[0]; i++) {
            for (int j = 0; j < s2[1]; j++) {
                m2[i][j] = rand() % 100;
            }
        }


        saveMatrixToFile(m1, m1_filename);
        saveMatrixToFile(m2, m2_filename);
        return 0;

    }

    else {

        std::vector<std::vector<int>> m1 = loadMatrixFromFile(m1_filename, s1[0], s1[1]);
        std::vector<std::vector<int>> m2 = loadMatrixFromFile(m2_filename, s2[0], s2[1]);

        clock_t start;
        start = clock();

        std::vector<std::vector<int>> result = matrixMultOpenMP(m1, m2);

        double time = double(clock() - start) / CLOCKS_PER_SEC;
        saveResult(result, time, s1[0], s1[1], s2[1]);

        //printMatrix(m1);
        //printMatrix(m2);
        //printMatrix(result);

        std::cout << std::endl << "Calculation time: " << time << " seconds." << std::endl;

        return 0;
    }
}
