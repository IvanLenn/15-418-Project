#include "linear_programming_seq.h"
#include <cassert>
#include <string>
#include <iostream>

LinearProgrammingSeq::LinearProgrammingSeq(const int n) : NumVar(n), NumCons(0) {}

void LinearProgrammingSeq::AddTarget(const std::vector<double>& T) {
    assert(T.size() == NumVar);
    Target.resize(NumVar);
    for (int i = 0; i < NumVar; i++) {
        Target[i] = T[i];
    }
}

void LinearProgrammingSeq::AddCons(const std::vector<std::vector<double>>& A) {
    for (auto &a : A) {
        assert(a.size() == NumVar + 1);
    }

    NumCons += A.size();
    for (auto &a : A) {
        Matrix.push_back(a);
    }
}

LinearProgrammingAnswer& LinearProgrammingSeq::Solve() const {}

void LinearProgrammingSeq::Print() const {
    auto f = [](const int i) {return "x_" + std::to_string(i);};
    std::cout << "Maximize ";
    for (int i = 0; i < NumVar; i++) {
        std::cout << Target[i] << "*" << f(i) << " ";
        if (i != NumVar - 1) {
            std::cout << "+ ";
        }
    }
    std::cout << '\n';

    std::cout << "Subject to: \n";
    for (int i = 0; i < NumCons; i++) {
        for (int j = 0; j < NumVar; j++) {
            std::cout << Matrix[i][j] << "*" << f(j) << " ";
            if (j != NumVar - 1) {
                std::cout << "+ ";
            }
        }
        std::cout << "<= " << Matrix[i][NumVar] << '\n';
    }
}