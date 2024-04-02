#include "linear_programming_seq.h"
#include <cassert>
#include <string>
#include <iostream>
#include <limits>

void LinearProgrammingAnswer::Print() const {
    switch (SolutionStatus) {
        case Infeasible:
            std::cout << "Infeasible\n";
            break;
        case Bounded:
            std::cout << "Bounded\n";
            std::cout << "Max: " << Max << '\n';
            std::cout << "Assignment: ";
            for (auto& a : Assignment) {
                std::cout << a << " ";
            }
            std::cout << '\n';
            break;
        case Unbounded:
            std::cout << "Unbounded\n";
            break;
    }
}

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

std::pair<int, int> LinearProgrammingSeq::FindPivot() {
    int pivot_row = -1;
    int pivot_col = -1;
    double min_ratio = std::numeric_limits<double>::max();
    double min = std::numeric_limits<double>::max();
    for (int i = 0; i < NumVar + NumCons; i++) {
        if (Tableau[NumCons][i] < min) {
            min = Tableau[NumCons][i];
            pivot_col = i;
        }
    }
    if (min >= 0) {
        // Answer = LinearProgrammingAnswer(Tableau[NumCons][NumVar + NumCons + 1],);
        Answer.SolutionStatus = LinearProgrammingAnswer::Status::Bounded;
        Answer.Max = Tableau[NumCons][NumVar + NumCons + 1];
        /*************************
        TODO: Assignment of variables not yet implemented
        *************************/
        return std::make_pair(pivot_row, pivot_col);
    }
    for (int i = 0; i < NumCons; i++) {
        if (Tableau[i][pivot_col] > 0) {
            double ratio = Tableau[i][NumVar + NumCons + 1] / Tableau[i][pivot_col];
            if (ratio < min_ratio) {
                min_ratio = ratio;
                pivot_row = i;
            }
        }
    }
    if (pivot_row == -1) {
        Answer.SolutionStatus = LinearProgrammingAnswer::Status::Unbounded;
        return std::make_pair(pivot_row, pivot_col);
    }
    return std::make_pair(pivot_row, pivot_col);
}

void LinearProgrammingSeq::Eliminate(int pivot_row, int pivot_col) {
    double pivot = Tableau[pivot_row][pivot_col];
    for (int i = 0; i < NumVar + NumCons + 2; i++) {
        Tableau[pivot_row][i] /= pivot;
    }
    for (int i = 0; i < NumCons + 1; i++) {
        if (i == pivot_row) {
            continue;
        }
        double ratio = Tableau[i][pivot_col];
        for (int j = 0; j < NumVar + NumCons + 2; j++) {
            Tableau[i][j] -= ratio * Tableau[pivot_row][j];
        }
    }
}

LinearProgrammingAnswer& LinearProgrammingSeq::Solve() {
    InitTableau();
    /*************************
    TODO: Case of no solution
    *************************/
    while (true) {
        auto [pivot_row, pivot_col] = FindPivot();
        if (pivot_row == -1) {
            break;
        }
        Eliminate(pivot_row, pivot_col);
    }
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