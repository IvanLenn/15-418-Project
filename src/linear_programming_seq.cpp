#include "linear_programming_seq.h"
#include <cassert>
#include <string>
#include <iostream>
#include <limits>

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

void LinearProgrammingSeq::PrintM(std::vector<std::vector<double>>& T) const {
    for (int i = 0; i < T.size(); i++) {
        for (int j = 0; j < T[i].size(); j++) {
            std::cout << T[i][j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

void LinearProgrammingSeq::InitTableau() {
    Tableau.clear();
    Tableau.resize(NumCons + 1);
    for (int i = 0; i < NumCons + 1; i++) {
        Tableau[i].resize(NumVar + NumCons + 2, 0);
    }

    for (int i = 0; i < NumCons; i++) {
        for (int j = 0; j < NumVar; j++) {
            Tableau[i][j] = Matrix[i][j];
        }
        Tableau[i][NumVar + i] = 1;
        Tableau[i][NumVar + NumCons + 1] = Matrix[i][NumVar];
    }

    for (int i = 0; i < NumVar; i++) {
        Tableau[NumCons][i] = -Target[i];
    }
    Tableau[NumCons][NumVar + NumCons] = 1;
}

std::pair<int, int> LinearProgrammingSeq::FindPivot() const {
    int pivot_row = -1;
    int pivot_col = -1;
    double min_ratio = std::numeric_limits<double>::max();
    int idx_max = -1;
    
    return std::make_pair(pivot_row, pivot_col);
}

LinearProgrammingAnswer& LinearProgrammingSeq::Solve() {
    InitTableau();
    return Answer;
}

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