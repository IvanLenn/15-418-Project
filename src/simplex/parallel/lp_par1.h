#pragma once
#include <vector>
#include <iostream>
#include <lp_par.h>

class LinearProgramming1 : public LinearProgramming {
private:
    const double EPI = 1e-6;
    int NumVar, NumCons;
    std::vector<std::vector<double>> Matrix;
    std::vector<double> Target;
    std::vector<std::vector<double>> Tableau;
    std::vector<int> Basic, NonBasic;
    LinearProgrammingAnswer Answer{};

    void Init() override;
    bool Feasible() override;
    std::pair<int, int> FindPivot() override;
    void Eliminate(const int pivot_row, const int pivot_col);
    void PrintM(const std::vector<std::vector<double>>& T) const;
public:
    LinearProgramming1() {};
    LinearProgramming1(const int n) {};
    void AddTarget(const std::vector<double>& T);
    void AddCons(const std::vector<std::vector<double>>& A);
    LinearProgrammingAnswer& Solve();
    void Check() const;
    void Print() const;
};