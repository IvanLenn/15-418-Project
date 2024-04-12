#include <vector>
#include <iostream>
#include <cassert>
#include "linear_programming_seq.h"

struct Edge {
    int from, to, cap;
    Edge(int from, int to, int cap) : from(from), to(to), cap(cap) {}
};

struct FlowAnswer {
    double Max;
    std::vector<double> Assignment;
    FlowAnswer(double Max, std::vector<double>& Assignment) : Max(Max), Assignment(Assignment) {}
};

class Flow {
    int n, m, s, t;
    std::vector<Edge> G;
    Flow(int n, int s, int t, std::vector<Edge> G) : n(n), s(s), t(t), G(G) {
        m = G.size();
        for (auto &e : G) {
            assert(e.cap >= 0);
            assert(e.from >= 0 && e.from < n);
            assert(e.to >= 0 && e.to < n);
            assert(e.from != e.to);
        }
    }

    FlowAnswer Solve() {
        std::vector<double> target(m, 0);
        std::vector<std::vector<double>> matrix(2 * (m + n - 2), std::vector<double>(m + 1, 0));

        std::vector<int> label(n);
        int start = 2 * m;
        for (int i = 0; i < n; i++) {
            if (i == s || i == t) continue;
            label[i] = start++;
        }

        for (int i = 0; i < m; i++) {
            if (G[i].from == s) target[i] = 1;
            if (G[i].to == s) target[i] = -1;
        }

        for (int i = 0; i < m; i++) {
            matrix[2 * i][i] = -1;
            matrix[2 * i + 1][i] = 1;
            matrix[2 * i + 1][m] = G[i].cap;
        }

        for (int i = 0; i < m; i++) {
            if (G[i].from != s && G[i].from != t) {
                matrix[label[G[i].from]][i] = 1;
                matrix[label[G[i].from] + 1][i] = -1;
            }
            if (G[i].to != s && G[i].to != t) {
                matrix[label[G[i].to]][i] = -1;
                matrix[label[G[i].to] + 1][i] = 1;
            }
        }

        LinearProgrammingSeq T(m);
        T.AddTarget(target);
        T.AddCons(matrix);
        auto ans = T.Solve();
        return FlowAnswer(ans.Max, ans.Assignment);
    }
};

int main() {
    return 0;
}