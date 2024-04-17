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
    LinearProgrammingAnswer(double Max, const std::vector<double>& Assignment, Status SolutionStatus) :
                                Max(Max), Assignment(Assignment), SolutionStatus(SolutionStatus) {};
    void Print() const;
    bool operator==(const LinearProgrammingAnswer& rhs) const;
};

class LinearProgramming {
private:
    const double EPI = 1e-6;
    int NumVar, NumCons;
    std::vector<std::vector<double>> Matrix;
    std::vector<double> Target;
    std::vector<std::vector<double>> Tableau;
    std::vector<int> Basic, NonBasic;
    LinearProgrammingAnswer Answer{};

    void Init();
    bool Feasible();
    std::pair<int, int> FindPivot();
    void Eliminate(const int pivot_row, const int pivot_col);
    void PrintM(const std::vector<std::vector<double>>& T) const;
public:
    LinearProgramming(const int n);
    void AddTarget(const std::vector<double>& T);
    void AddCons(const std::vector<std::vector<double>>& A);
    LinearProgrammingAnswer& Solve();
    void Check() const;
    void Print() const;
};