#include <iostream>
#include <vector>
#include "linear_programming_seq.h"

int main() {
    LinearProgrammingSeq lps(2);
    std::vector<double> target = {1, 1};
    std::vector<std::vector<double>> cons = {{1, 1, 2}, {1, 2, 4}};
    lps.AddTarget(target);
    lps.AddCons(cons);
    lps.Print();
}