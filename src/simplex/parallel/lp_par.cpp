#include "lp_par.h"
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

bool LinearProgrammingAnswer::operator==(const LinearProgrammingAnswer& rhs) const {
    if (SolutionStatus != rhs.SolutionStatus) {
        return false;
    }
    if (SolutionStatus == Infeasible) {
        return true;
    }
    if (SolutionStatus == Unbounded) {
        return true;
    }
    if (std::abs(Max - rhs.Max) > 1e-6) {
        return false;
    }
    return true;
}