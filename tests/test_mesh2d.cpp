#include <gtest/gtest.h>
#include "../src/mesh2d.hpp"
#include <vector>
#include <cmath>

// Function to compute the error between u_exact and u_h
double computeError(const std::vector<double>& u_exact, const std::vector<double>& u_h) {
    double error = 0.0;
    for (size_t i = 0; i < u_exact.size(); ++i) {
        error += std::pow(u_exact[i] - u_h[i], 2);
    }
    return std::sqrt(error);
}

TEST(Mesh2DTest, ComputeError) {
    std::vector<double> u_exact = {1.0, 2.0, 3.0}; // Example exact values
    std::vector<double> u_h = {0.9, 2.1, 2.9}; // Example computed values

    double error = computeError(u_exact, u_h);
    EXPECT_NEAR(error, 0.173, 1e-3); // Example tolerance for the error
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
