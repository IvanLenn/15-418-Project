#include <cassert>
#include "flow_par.h"
#include <lp_par1.h>
#include <mpi.h>

/*************************
PRECONDITION: Only process nproc - 1 inputs the useful Graph. All others are ignored
*************************/
Flow::Flow(int n, int s, int t, std::vector<Edge> G) : n(n), s(s), t(t), G(G) {
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    m = G.size();
    for (auto &e : G) {
        assert(e.cap >= 0);
        assert(e.from >= 0 && e.from < n);
        assert(e.to >= 0 && e.to < n);
    }
}

FlowAnswer Flow::Solve() {
    MPI_Bcast(&m, 1, MPI_INT, nproc - 1, MPI_COMM_WORLD);
    LinearProgramming1 T(m);

    if (pid ==  nproc - 1) {
        std::vector<double> target(m, 0);
        std::vector<std::vector<double>> matrix(2 * (m + n - 2), std::vector<double>(m + 1, 0));

        std::vector<int> label(n);
        int start = 2 * m;
        for (int i = 0; i < n; i++) {
            if (i == s || i == t) continue;
            label[i] = start;
            start += 2;
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

        T.AddTarget(target);
        T.AddCons(matrix);
    }
    auto ans = *T.Solve();
    // return FlowAnswer(ans.Max, ans.Assignment);
    auto t = std::vector<double>(1, 1);
    return FlowAnswer(1, t);
}