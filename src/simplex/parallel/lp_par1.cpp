#include "lp_par1.h"
#include <cassert>
#include <string>
#include <iostream>
#include <limits>
#include <mpi.h>

LinearProgramming1::LinearProgramming1() {
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    NumVar = 0;
    NumCons = 0;
    Target = nullptr;
    MatrixData = nullptr;
    Matrix = nullptr;
    TableauData = nullptr;
    Tableau = nullptr;
    Basic = nullptr;
    NonBasic = nullptr;
}

LinearProgramming1::LinearProgramming1(const int n) {
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    NumVar = n;
    NumCons = 0;
    Target = nullptr;
    MatrixData = nullptr;
    Matrix = nullptr;
    TableauData = nullptr;
    Tableau = nullptr;
    Basic = nullptr;
    NonBasic = nullptr;
}

LinearProgramming1::~LinearProgramming1() {
    delete[] Target;
    delete[] Matrix;
    delete[] MatrixData;
    delete[] Tableau;
    delete[] TableauData;
    delete[] Basic;
    delete[] NonBasic;
}

/*************************
PRECONDITION: Only pid nproc - 1 will call this function
*************************/
void LinearProgramming1::AddTarget(const std::vector<double>& T) {
    assert(pid == nproc - 1);
    assert(T.size() == NumVar);
    Target = new double[NumVar];
    for (int i = 0; i < NumVar; i++) {
        Target[i] = T[i];
    }
}

/*************************
PRECONDITION: Only pid nproc - 1 will call this function
*************************/
void LinearProgramming1::AddCons(const std::vector<std::vector<double>>& A) {
    assert(pid == nproc - 1);
    for (auto &a : A) {
        assert(a.size() == NumVar + 1);
    }

    NumCons += A.size();
    int pad = (NumCons + nproc - 1) / nproc * nproc - NumCons;
    double* bData = new double[(NumCons + pad) * (NumVar + 1)];
    double** b = new double*[NumCons + pad];
    for (int i = 0; i < NumCons + pad; i++) {
        b[i] = bData + i * (NumVar + 1);
    }
    for (int i = 0; i < NumCons - A.size(); i++) {
        for (int j = 0; j <= NumVar; j++) {
            b[i][j] = Matrix[i][j];
        }
    }
    for (int i = NumCons - A.size(); i < NumCons; i++) {
        for (int j = 0; j <= NumVar; j++) {
            b[i][j] = A[i - NumCons + A.size()][j];
        }
    }
    for (int i = NumCons; i < NumCons + pad; i++) {
        for (int j = 0; j <= NumVar; j++) {
            b[i][j] = 0.0f;
        }
    }
    delete[] MatrixData;
    delete[] Matrix;
    MatrixData = bData;
    Matrix = b;
}

/*************************
PRECONDITION: Only calling the whole function once; Also once for Solve;
*************************/
void LinearProgramming1::Init() {
    // Broadcast the number of constraints
    if (pid == nproc - 1) {
        MPI_Bcast(&NumCons, 1, MPI_INT, nproc - 1, MPI_COMM_WORLD);
    }

    // Allocate all required arrays across processes
    int task_required = (NumCons + nproc - 1) / nproc;
    StartCons = pid * task_required;
    EndCons = std::min(NumCons, (pid + 1) * task_required);
    TableauData = new double[task_required * (NumVar + 1)];
    Tableau = new double*[task_required];
    for (int j = 0; j < task_required; j++) {
        Tableau[j] = TableauData + j * (NumVar + 1);
    }
    Basic = new int[NumCons];
    NonBasic = new int[NumVar];
    for (int i = 0; i < NumCons; i++) {
        Basic[i] = NumVar + i;
    }
    for (int i = 0; i < NumVar; i++) {
        NonBasic[i] = i;
    }
    Answer = LinearProgrammingAnswer();
    MPI_Scatter(MatrixData, task_required * (NumVar + 1), MPI_DOUBLE, TableauData, task_required * (NumVar + 1), MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);

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
    auto t = LinearProgrammingAnswer();
    return t;
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
    if (pid == 0) std::cout << "Maximize ";
    if (pid == 0) {
        for (int i = 0; i < NumVar; i++) {
            std::cout << Target[i] << "*" << f(i) << " ";
            if (i != NumVar - 1) {
                std::cout << "+ ";
            }
        }
        std::cout << '\n';
    }

    if (pid == 0) std::cout << "Subject to: \n";
    for (int i = 0; i < nproc; i++) {
        if (pid == i) {
            for (int j = StartCons; j < EndCons; j++) {
                for (int k = 0; k < NumVar; k++) {
                    std::cout << Matrix[j][k] << "*" << f(k) << " ";
                    if (k != NumVar - 1) {
                        std::cout << "+ ";
                    }
                }
                std::cout << "<= " << Matrix[j][NumVar] << '\n';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    std::cout << '\n';
}