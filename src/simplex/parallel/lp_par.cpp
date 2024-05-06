#include "lp_par.h"
#include <cassert>
#include <string>
#include <iostream>
#include <limits>
#include <mpi.h>

LinearProgramming::LinearProgramming(const int n) {
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
    buffer = nullptr;
    Answer = new LinearProgrammingAnswer();
}

LinearProgramming::~LinearProgramming() {
    delete[] Target;
    delete[] Matrix;
    delete[] MatrixData;
    delete[] Tableau;
    delete[] TableauData;
    delete[] Basic;
    delete[] NonBasic;
    delete[] buffer;
    delete Answer;
}

void LinearProgramming::AddTarget(const std::vector<double>& T) {
    assert(pid == nproc - 1);
    assert(T.size() == NumVar);
    InputTarget.resize(NumVar);
    for (int i = 0; i < NumVar; i++) {
        InputTarget[i] = T[i];
    }
    Target = new double[NumVar + 1];
    for (int i = 0; i < NumVar; i++) {
        Target[i] = T[i];
    }
    Target[NumVar] = 0;
}

void LinearProgramming::AddCons(const std::vector<std::vector<double>>& A) {
    assert(pid == nproc - 1);
    for (auto &a : A) {
        assert(a.size() == NumVar + 1);
    }
    for (auto a : A) {
        InputMatrix.push_back(a);
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

void LinearProgramming::Init() {
    // Broadcast the number of constraints
    MPI_Bcast(&NumCons, 1, MPI_INT, nproc - 1, MPI_COMM_WORLD);
    if (pid != nproc - 1) {
        Target = new double[NumVar + 1];
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
    buffer = new double[std::max(NumVar + 1, task_required * nproc)];
    MPI_Scatter(MatrixData, task_required * (NumVar + 1), MPI_DOUBLE, TableauData, task_required * (NumVar + 1), MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);
    MPI_Bcast(Target, NumVar + 1, MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);
}

std::pair<int, int> LinearProgramming::FindPivot() {
    int pivot_row = -1;
    int pivot_col = -1;
    double p = 0.0f;
    for (int i = 0; i < NumVar; i++) {
        if (Target[i] > p) {
            pivot_col = i;
            p = Target[i];
        }
    }

    if (p < EPI) {
        int task_required = (NumCons + nproc - 1) / nproc;
        double* tmp = new double[task_required];
        for (int i = 0; i < task_required; i++) {
            tmp[i] = Tableau[i][NumVar];
        }
        MPI_Gather(tmp, task_required, MPI_DOUBLE, buffer, task_required, MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);
        if (pid == nproc - 1) {
            Answer -> SolutionStatus = LinearProgrammingAnswer::Bounded;
            Answer -> Max = -Target[NumVar];
            Answer -> Assignment.resize(NumVar, 0);
            for (int i = 0; i < NumVar; i++) {
                if (Basic[i] < NumVar) {
                    Answer -> Assignment[Basic[i]] = buffer[i];
                }
            }
        }
        delete[] tmp;
        return std::make_pair(pivot_row, pivot_col);
    }

    p = std::numeric_limits<double>::max();
    bool unbound = true;
    for (int i = StartCons; i < EndCons; i++) {
        if (Tableau[i - StartCons][pivot_col] > EPI) {
            double ratio = Tableau[i - StartCons][NumVar] / Tableau[i - StartCons][pivot_col];
            if (ratio < p) {
                unbound = false;
                pivot_row = i;
                p = ratio;
            }
        }
    }

    in.r = p;
    in.idx = pivot_row;
    MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    if (out.idx >= 0) unbound = false; 
    pivot_row = out.idx;

    if (unbound) {
        Answer -> SolutionStatus = LinearProgrammingAnswer::Unbounded;
        return std::make_pair(pivot_row, pivot_col);
    }
    return std::make_pair(pivot_row, pivot_col);
}

void LinearProgramming::Eliminate(const int pivot_row, const int pivot_col) {
    std::swap(Basic[pivot_row], NonBasic[pivot_col]);
    int task_required = (NumCons + nproc - 1) / nproc;
    int pid_row_local = pivot_row / task_required;
    if (StartCons <= pivot_row && pivot_row < EndCons) {
        Tableau[pivot_row - StartCons][pivot_col] = 1.0 / Tableau[pivot_row - StartCons][pivot_col];
        for (int i = 0; i <= NumVar; i++) {
            if (i != pivot_col) {
                Tableau[pivot_row - StartCons][i] *= Tableau[pivot_row - StartCons][pivot_col];
            }
        }
    }
    if (pid == pid_row_local) {
        for (int i = 0; i <= NumVar; i++) {
            buffer[i] = Tableau[pivot_row - StartCons][i];
        }
    }
    MPI_Bcast(buffer, NumVar + 1, MPI_DOUBLE, pid_row_local, MPI_COMM_WORLD);

    for (int i = StartCons; i < EndCons; i++) {
        if (i != pivot_row) {
            for (int j = 0; j <= NumVar; j++) {
                if (j != pivot_col) {
                    Tableau[i - StartCons][j] -= Tableau[i - StartCons][pivot_col] * buffer[j];
                }
            }
            Tableau[i - StartCons][pivot_col] = -Tableau[i - StartCons][pivot_col] * buffer[pivot_col];
        }
    }
    for (int j = 0; j <= NumVar; j++) {
        if (j != pivot_col) {
            Target[j] -= Target[pivot_col] * buffer[j];
        }
    }
    Target[pivot_col] = -Target[pivot_col] * buffer[pivot_col];
}


bool LinearProgramming::Feasible() {
    int pivot_row = -1;
    int pivot_col = -1;
    while (true) {
        double p = std::numeric_limits<double>::max();
        for (int i = StartCons; i < EndCons; i++) {
            if (Tableau[i - StartCons][NumVar] < p) {
                pivot_row = i;
                p = Tableau[i - StartCons][NumVar];
            }
        }
        in.r = p;
        in.idx = pivot_row;
        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

        if (out.r > -EPI) {
            return true;
        }
        pivot_row = out.idx;

        p = 0.0f;
        int pid_row_local = pivot_row / (EndCons - StartCons);
        if (pid == pid_row_local) {
            for (int i = 0; i < NumVar; i++) {
                if (Tableau[pivot_row - StartCons][i] < p) {
                    pivot_col = i;
                    p = Tableau[pivot_row - StartCons][i];
                }
            }
        }
        MPI_Bcast(&pivot_col, 1, MPI_INT, pid_row_local, MPI_COMM_WORLD);
        if (pivot_col == -1) {
            return false;
        }

        p = std::numeric_limits<double>::max();
        if (pid == pid_row_local) {
            for (int i = pivot_row - StartCons; i < EndCons - StartCons; i++) {
                if (Tableau[i][pivot_col] > EPI) {
                    double tmp = Tableau[i][NumVar] / Tableau[i][pivot_col];
                    if (tmp < p) {
                        pivot_row = i + StartCons;
                        p = tmp;
                    }
                }
            }
        }
        else if (pid > pid_row_local) {
            for (int i = StartCons; i < EndCons; i++) {
                if (Tableau[i - StartCons][pivot_col] > EPI) {
                    double tmp = Tableau[i - StartCons][NumVar] / Tableau[i - StartCons][pivot_col];
                    if (tmp < p) {
                        pivot_row = i;
                        p = tmp;
                    }
                }
            }
        }
        in.r = p;
        in.idx = pivot_row;
        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        Eliminate(pivot_row, pivot_col);
    }
}

LinearProgrammingAnswer* LinearProgramming::Solve() {
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

void LinearProgramming::Check() const {
    // Broadcast Answer
    MPI_Bcast(&Answer -> SolutionStatus, 1, MPI_INT, nproc - 1, MPI_COMM_WORLD);
    MPI_Bcast(&Answer -> Max, 1, MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);
    MPI_Bcast(Answer -> Assignment.data(), NumVar, MPI_DOUBLE, nproc - 1, MPI_COMM_WORLD);

    if (Answer -> SolutionStatus != LinearProgrammingAnswer::Bounded) {
        return;
    }

    for (int i = StartCons; i < EndCons; i++) {
        double sum = 0.0f;
        for (int j = 0; j < NumVar; j++) {
            sum += Answer -> Assignment[j] * Tableau[i - StartCons][j];
        }
        if (sum > Tableau[i - StartCons][NumVar] + EPI) {
            std::cerr << "Check failed on " << i << " th constraint\n";
            exit(EXIT_FAILURE);
        }
    }

    if (pid == nproc - 1) {
        double max = 0.0f;
        for (int i = 0; i < NumVar; i++) {
            max += Answer -> Assignment[i] * Target[i];
        }
        if (max < Answer -> Max - EPI || max > Answer -> Max + EPI) {
            std::cerr << "Check failed\n";
            exit(EXIT_FAILURE);
        }
    }
}

void LinearProgramming::Print() const {
    if (pid != nproc - 1) return;
    auto f = [](const int i) {return "x_" + std::to_string(i);};
    std::cout << "Maximize ";
    for (int i = 0; i < NumVar; i++) {
        std::cout << InputTarget[i] << "*" << f(i) << " ";
        if (i != NumVar - 1) {
            std::cout << "+ ";
        }
    }
    std::cout << '\n';
    std::cout << "Subject to: \n";
    for (int i = 0; i < NumCons; i++) {
        for (int j = 0; j < NumVar; j++) {
            std::cout << InputMatrix[i][j] << "*" << f(j) << " ";
            if (j != NumVar - 1) {
                std::cout << "+ ";
            }
        }
        std::cout << "<= " << InputMatrix[i][NumVar] << '\n';
    }
}

void LinearProgramming::dbg() const {
    std::cout << "From pid: " << pid << " in dbg: " << StartCons << ' ' << EndCons << '\n';
    if (pid == nproc - 1) std::cout << "Tableau:\n";
    for (int i = 0; i < nproc; i++) {
        if (pid == i) {
            for (int j = 0; j < EndCons - StartCons; j++) {
                for (int k = 0; k <= NumVar; k++) {
                    std::cout << Tableau[j][k] << " ";
                }
                std::cout << '\n';
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (pid == nproc - 1) {
        for (int i = 0; i <= NumVar; i++) {
            std::cout << Target[i] << " ";
        }
        std::cout << '\n';
        std::cout << "Basic: ";
        for (int i = 0; i < NumCons; i++) {
            std::cout << Basic[i] << " ";
        }
        std::cout << '\n';
        std::cout << "NonBasic: ";
        for (int i = 0; i < NumVar; i++) {
            std::cout << NonBasic[i] << " ";
        }
        std::cout << '\n';
        Answer -> Print();
        std::cout << '\n';
    }
}