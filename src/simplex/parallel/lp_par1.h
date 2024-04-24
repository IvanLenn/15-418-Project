#pragma once
#include <vector>
#include <iostream>

struct LinearProgrammingAnswer {
    enum Status : int {
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
    std::vector<std::vector<double>> InputMatrix;
    std::vector<double> InputTarget;
    double* MatrixData;
    double** Matrix;
    double* Target;
    double* TableauData;
    double** Tableau;
    double* buffer;
    int* Basic, *NonBasic;
    LinearProgrammingAnswer* Answer;
    struct {
        double r;
        int idx;
    } in, out;

    void Init();
    bool Feasible();
    std::pair<int, int> FindPivot();
    void Eliminate(const int pivot_row, const int pivot_col);
public:
    /**
     * Destructor for LinearProgramming1.
     * Cleans up resources used by the instance.
     */
    ~LinearProgramming1();

    /**
     * Constructor for LinearProgramming1.
     * Initializes a LinearProgramming1 instance with specified variable capacity.
     * @param n The number of variables.
     */
    LinearProgramming1(const int n);
    
    /**
     * Adds the target function coefficients to the linear program.
     * @param T A vector containing the coefficients of the target function.
     * PRECONDITION: Only processor nproc - 1 can call this function; Others raise an exception.
     * PRECONDITION: T.size() == NumVar.
     */
    void AddTarget(const std::vector<double>& T);

    /**
     * Adds the constraints coefficients to the linear program.
     * @param A A 2-Dimensional vector containing the coefficients of constraints.
     * PRECONDITION: Only processor nproc - 1 can call this function; Others raise an exception.
     * PRECONDITION: Any A[i].size() == NumVar + 1.
     */
    void AddCons(const std::vector<std::vector<double>>& A);

    /**
     * Solves the linear program.
     * @return A LinearProgrammingAnswer object containing the solution.
     * PRECONDITION: Can only be called once. [TO DO]
     * POSTCONDITION: Only return from processor nproc - 1 is valid.
     */
    LinearProgrammingAnswer* Solve();

    /**
     * Checks the validity of the LinearProgramming1 instance.
     * Raises an exception if the instance is invalid.
     */
    void Check() const;

    /**
     * Prints the LinearProgramming1 instance.
     */
    void Print() const;

    /**
     * Debugging function to print the LinearProgramming1 instance.
     */
    void dbg() const;
};