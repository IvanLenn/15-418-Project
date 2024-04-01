#pragma once
#include <vector>

struct LinearProgrammingAnswer {
    double Max;
    std::vector<double> Assignment;
    LinearProgrammingAnswer(const double Max, std::vector<double>& Assignment) :
                                Max(Max), Assignment(Assignment) {};
};

class LinearProgrammingSeq {
private:
    int NumVar, NumCons;
    std::vector<std::vector<double>> Matrix;
    std::vector<double> Target;
public:
    LinearProgrammingSeq(const int n);
    void AddTarget(const std::vector<double>& T);
    void AddCons(const std::vector<std::vector<double>>& A);
    LinearProgrammingAnswer& Solve() const;
    void Print() const;
};