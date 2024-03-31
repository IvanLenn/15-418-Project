#pragma once
#include <vector>

struct LinearProgrammingAnswer {
    double Max;
    std::vector<double> Assignment;
    void LinearProgrammingAnswer(const double Max, std::vector<double>& Assignment) :
                                Max(Max), Assignment(Assignment);
}

class LinearProgrammingSeq {
private:
    int NumVar, NumCons;
    std::vector<std::vector<double>> Matrix;
public:
    void LinearProgrammingSeq(const int n);
    void AddCons(const int m, const std::vector<std::vector<double>>& A,
                 const std::vector<std::vector<double>>& Target);
    LinearProgrammingAnswer& Solve() const;
}