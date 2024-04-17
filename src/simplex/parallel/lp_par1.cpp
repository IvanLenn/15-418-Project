#include "lp_par1.h"
#include <cassert>
#include <string>
#include <iostream>
#include <limits>
#include <omp.h>

void LinearProgramming1::AddTarget(const std::vector<double>& T) {
    assert(T.size() == NumVar);
    Target.resize(NumVar);
    for (int i = 0; i < NumVar; i++) {
        Target[i] = T[i];
    }
}

void LinearProgramming1::AddCons(const std::vector<std::vector<double>>& A) {
    for (auto &a : A) {
        assert(a.size() == NumVar + 1);
    }

    NumCons += A.size();
    for (auto &a : A) {
        Matrix.push_back(a);
    }
}

void LinearProgramming1::PrintM(const std::vector<std::vector<double>>& T) const {
    for (int i = 0; i < T.size(); i++) {
        for (int j = 0; j < T[i].size(); j++) {
            std::cout << T[i][j] << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

void LinearProgramming1::Init() {
    Tableau.clear();
    Basic.clear();
    NonBasic.clear();
    Answer = LinearProgrammingAnswer();
    Tableau.resize(NumCons + 1, std::vector<double>(NumVar + 1, 0));

    for (int i = 0; i < NumCons; i++) {
        for (int j = 0; j <= NumVar; j++) {
            Tableau[i][j] = Matrix[i][j];
        }
    }

    for (int i = 0; i < NumVar; i++) {
        Tableau[NumCons][i] = Target[i];
    }

    for (int i = 0; i < NumCons; i++) {
        Basic.push_back(NumVar + i);
    }
    for (int i = 0; i < NumVar; i++) {
        NonBasic.push_back(i);
    }
}

std::pair<int, int> LinearProgramming1::FindPivot() {
    int pivot_row = -1;
    int pivot_col = -1;
    double p = 0.0f;
    for (int i = 0; i < NumVar; i++) {
        if (Tableau[NumCons][i] > p) {
            pivot_col = i;
            p = Tableau[NumCons][i];
        }
    }

    if (p < EPI) {
        Answer.SolutionStatus = LinearProgrammingAnswer::Bounded;
        Answer.Max = -Tableau[NumCons][NumVar];
        Answer.Assignment.resize(NumVar, 0);
        for (int i = 0; i < NumVar; i++) {
            if (Basic[i] < NumVar) {
                Answer.Assignment[Basic[i]] = Tableau[i][NumVar];
            }
        }
        return std::make_pair(pivot_row, pivot_col);
    }

    p = std::numeric_limits<double>::max();
    bool unbound = true;
    for (int i = 0; i < NumCons; i++) {
        if (Tableau[i][pivot_col] > EPI) {
            double ratio = Tableau[i][NumVar] / Tableau[i][pivot_col];
            if (ratio < p) {
                unbound = false;
                pivot_row = i;
                p = ratio;
            }
        }
    }
    if (unbound) {
        Answer.SolutionStatus = LinearProgrammingAnswer::Unbounded;
        return std::make_pair(pivot_row, pivot_col);
    }
    return std::make_pair(pivot_row, pivot_col);
}

void LinearProgramming1::Eliminate(const int pivot_row, const int pivot_col) {
    std::swap(Basic[pivot_row], NonBasic[pivot_col]);
    Tableau[pivot_row][pivot_col] = 1.0 / Tableau[pivot_row][pivot_col];
    for (int i = 0; i <= NumVar; i++) {
        if (i != pivot_col) {
            Tableau[pivot_row][i] *= Tableau[pivot_row][pivot_col];
        }
    }
    for (int i = 0; i <= NumCons; i++) {
        if (i != pivot_row) {
            for (int j = 0; j <= NumVar; j++) {
                if (j != pivot_col) {
                    Tableau[i][j] -= Tableau[i][pivot_col] * Tableau[pivot_row][j];
                }
            }
            Tableau[i][pivot_col] = -Tableau[i][pivot_col] * Tableau[pivot_row][pivot_col];
        }
    }
}


bool LinearProgramming1::Feasible() {
    int pivot_row = -1;
    int pivot_col = -1;
    while (true) {
        double p = std::numeric_limits<double>::max();
        for (int i = 0; i < NumCons; i++) {
            if (Tableau[i][NumVar] < p) {
                pivot_row = i;
                p = Tableau[i][NumVar];
            }
        }
        if (p > -EPI) {
            return true;
        }
        p = 0.0f;
        for (int i = 0; i < NumVar; i++) {
            if (Tableau[pivot_row][i] < p) {
                pivot_col = i;
                p = Tableau[pivot_row][i];
            }
        }
        if (p > -EPI) {
            return false;
        }
        p = Tableau[pivot_row][NumVar] / Tableau[pivot_row][pivot_col];
        for (int i = pivot_row + 1; i < NumCons; i++) {
            if (Tableau[i][pivot_col] > EPI) {
                double tmp = Tableau[i][NumVar] / Tableau[i][pivot_col];
                if (tmp < p) {
                    pivot_row = i;
                    p = tmp;
                }
            }
        }
        Eliminate(pivot_row, pivot_col);
    }
}

LinearProgrammingAnswer& LinearProgramming1::Solve() {
    Init();
    if (!Feasible()) {
        return Answer;
    }
    while (true) {
        auto [pivot_row, pivot_col] = FindPivot();
        if (pivot_row == -1) {
            break;
        }
        Eliminate(pivot_row, pivot_col);
    }
    return Answer;
}

void LinearProgramming1::Check() const {
    if (Answer.SolutionStatus != LinearProgrammingAnswer::Bounded) {
        return;
    }
    double max = 0.0f;
    for (int i = 0; i < NumCons; i++) {
        double sum = 0.0f;
        for (int j = 0; j < NumVar; j++) {
            sum += Answer.Assignment[j] * Matrix[i][j];
        }
        if (sum > Matrix[i][NumVar] + EPI) {
            std::cerr << "Check failed\n";
            exit(EXIT_FAILURE);
        }
    }
    for (int i = 0; i < NumVar; i++) {
        max += Answer.Assignment[i] * Target[i];
    }
    if (max < Answer.Max - EPI) {
        std::cerr << "Check failed\n";
        exit(EXIT_FAILURE);
    }
}

void LinearProgramming1::Print() const {
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