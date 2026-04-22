#ifndef SRC_HPP
#define SRC_HPP

#include <iostream>
#include <string>
#include <exception>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <utility>
#include "fraction.hpp"

// 如果你不需要使用 matrix 类，请将 IGNORE_MATRIX 改为 0
// #define IGNORE_MATRIX 0
#define IGNORE_MATRIX 0

#if IGNORE_MATRIX

class matrix {
private:

    int m, n;
    fraction **data;

public:

    matrix() {
        m = n = 0;
        data = nullptr;
    }

    matrix(int m_, int n_);
    matrix(const matrix &obj);
    matrix(matrix &&obj) noexcept;
    ~matrix();
    matrix &operator=(const matrix &obj);
    fraction &operator()(int i, int j);
    friend matrix operator*(const matrix &lhs, const matrix &rhs);
    matrix transposition();
    fraction determination();
};

#endif

class resistive_network {
private:
    int n, m;
    struct Edge { int u, v; fraction g; };
    std::vector<Edge> edges;
    std::vector<std::vector<fraction>> M; // (n-1)x(n-1) reduced Laplacian

    static std::vector<fraction> solve_linear(std::vector<std::vector<fraction>> A,
                                              const std::vector<fraction> &b) {
        int N = (int)A.size();
        if (N == 0) throw resistive_network_error();
        std::vector<std::vector<fraction>> aug(N, std::vector<fraction>(N + 1));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) aug[i][j] = A[i][j];
            aug[i][N] = b[i];
        }
        for (int col = 0, row = 0; col < N && row < N; ++col, ++row) {
            int piv = row;
            for (int i = row; i < N; ++i) {
                if (!(aug[i][col] == fraction(0))) { piv = i; break; }
            }
            if (aug[piv][col] == fraction(0)) {
                throw resistive_network_error();
            }
            if (piv != row) std::swap(aug[piv], aug[row]);
            fraction pivval = aug[row][col];
            for (int j = col; j <= N; ++j) aug[row][j] = aug[row][j] / pivval;
            for (int i = 0; i < N; ++i) if (i != row) {
                fraction factor = aug[i][col];
                if (factor == fraction(0)) continue;
                for (int j = col; j <= N; ++j) {
                    aug[i][j] = aug[i][j] - factor * aug[row][j];
                }
            }
        }
        std::vector<fraction> x(N);
        for (int i = 0; i < N; ++i) x[i] = aug[i][N];
        return x;
    }

public:

    resistive_network(int interface_size_, int connection_size_, int from[], int to[], fraction resistance[]) {
        n = interface_size_;
        m = connection_size_;
        edges.reserve(m);
        std::vector<std::vector<fraction>> L(n, std::vector<fraction>(n, fraction(0)));
        for (int k = 0; k < m; ++k) {
            int u = from[k];
            int v = to[k];
            fraction g = fraction(1) / resistance[k];
            edges.push_back({u, v, g});
            L[u-1][u-1] = L[u-1][u-1] + g;
            L[v-1][v-1] = L[v-1][v-1] + g;
            L[u-1][v-1] = L[u-1][v-1] - g;
            L[v-1][u-1] = L[v-1][u-1] - g;
        }
        if (n >= 2) {
            M.assign(n-1, std::vector<fraction>(n-1, fraction(0)));
            for (int i = 0; i < n-1; ++i)
                for (int j = 0; j < n-1; ++j)
                    M[i][j] = L[i][j];
        } else {
            M.clear();
        }
    }

    ~resistive_network() = default;

    fraction get_equivalent_resistance(int interface_id1, int interface_id2) {
        int a = interface_id1;
        int b = interface_id2;
        if (n <= 1) throw resistive_network_error();
        std::vector<fraction> bvec(n-1, fraction(0));
        if (a < n) bvec[a-1] = bvec[a-1] + fraction(1);
        if (b < n) bvec[b-1] = bvec[b-1] - fraction(1);
        std::vector<fraction> x = solve_linear(M, bvec);
        if (a < n && b < n) return x[a-1] - x[b-1];
        if (a < n && b == n) return x[a-1];
        if (a == n && b < n) return fraction(0) - x[b-1];
        return fraction(0);
    }

    fraction get_voltage(int id, fraction current[]) {
        if (id >= n || id < 1) throw resistive_network_error();
        if (n == 1) return fraction(0);
        std::vector<fraction> bvec(n-1, fraction(0));
        for (int i = 0; i < n-1; ++i) bvec[i] = current[i];
        std::vector<fraction> x = solve_linear(M, bvec);
        return x[id-1];
    }

    fraction get_power(fraction voltage[]) {
        fraction total(0);
        for (const auto &e : edges) {
            fraction du = voltage[e.u - 1] - voltage[e.v - 1];
            total = total + e.g * du * du;
        }
        return total;
    }
};


#endif //SRC_HPP
