#pragma once
#include <vector>

struct LinearProgrammingAnswer {
    enum Status {
        Infeasible,
        Bounded,
        Unbounded
    };
    double Max;
    std::vector<double> Assignment;
    Status SolutionStatus;
    LinearProgrammingAnswer() : Max(0), Assignment(std::vector<double>()), SolutionStatus(Infeasible) {};
    LinearProgrammingAnswer(double Max, std::vector<double>& Assignment, Status SolutionStatus) :
                                Max(Max), Assignment(Assignment), SolutionStatus(SolutionStatus) {};
    void Print() const;
};

class LinearProgrammingSeq {
private:
    int NumVar, NumCons;
    std::vector<std::vector<double>> Matrix;
    std::vector<double> Target;
    std::vector<std::vector<double>> Tableau;
    LinearProgrammingAnswer Answer{};

    void InitTableau();
    std::pair<int, int> FindPivot();
    void Eliminate(int pivot_row, int pivot_col);
    void PrintM(std::vector<std::vector<double>>& T) const;
public:
    LinearProgrammingSeq(const int n);
    void AddTarget(const std::vector<double>& T);
    void AddCons(const std::vector<std::vector<double>>& A);
    LinearProgrammingAnswer& Solve();
    void Print() const;
};