#include <iostream>
#include <vector>
#include "linear_programming_seq.h"

int main() {
    LinearProgrammingSeq lps(3);
    std::vector<double> target = {4, 1, 4};
    std::vector<std::vector<double>> cons = {{2, 1, 1, 2}, {1, 2, 3, 4}, {2, 2, 1, 8}};
    lps.AddTarget(target);
    lps.AddCons(cons);
    lps.Print();
    LinearProgrammingAnswer Ans = lps.Solve();
    Ans.Print();
}