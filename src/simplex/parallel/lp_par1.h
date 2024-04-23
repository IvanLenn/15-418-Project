#pragma once
#include <vector>
#include <iostream>

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
    void Print() const {
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

    bool operator!=(const LinearProgrammingAnswer& rhs) const {
        if (SolutionStatus != rhs.SolutionStatus) {
            return true;
        }
        if (SolutionStatus == Infeasible || SolutionStatus == Unbounded) {
            return false;
        }
        if (std::abs(Max - rhs.Max) > 1e-6) {
            return true;
        }
        return false;
    }
};

class LinearProgramming1 {
private:
    int pid, nproc;
    int StartCons, EndCons;

    const double EPI = 1e-6;
    int NumVar, NumCons;
    double* MatrixData;
    double** Matrix;
    double* Target;
    double* TableauData;
    double** Tableau;
    int* Basic, *NonBasic;
    LinearProgrammingAnswer Answer{};

    void Init();
    bool Feasible();
    std::pair<int, int> FindPivot();
    void Eliminate(const int pivot_row, const int pivot_col);
public:
    LinearProgramming1();
    ~LinearProgramming1();
    LinearProgramming1(const int n);
    void AddTarget(const std::vector<double>& T);
    void AddCons(const std::vector<std::vector<double>>& A);
    LinearProgrammingAnswer& Solve();
    void Check() const;
    void Print() const;
};