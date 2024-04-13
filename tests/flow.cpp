#include <vector>
#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "linear_programming_seq.h"

struct Edge {
    int from, to;
    double cap;
    Edge(int from, int to, double cap) : from(from), to(to), cap(cap) {}
};

struct FlowAnswer {
    double Max;
    std::vector<double> Assignment;
    FlowAnswer(double Max, std::vector<double>& Assignment) : Max(Max), Assignment(Assignment) {}
};

class Flow {
    int n, m, s, t;
    std::vector<Edge> G;
public:
    Flow(int n, int s, int t, std::vector<Edge> G) : n(n), s(s), t(t), G(G) {
        m = G.size();
        for (auto &e : G) {
            assert(e.cap >= 0);
            assert(e.from >= 0 && e.from < n);
            assert(e.to >= 0 && e.to < n);
        }
    }

    FlowAnswer Solve() {
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

        LinearProgrammingSeq T(m);
        T.AddTarget(target);
        T.AddCons(matrix);
        auto ans = T.Solve();
        return FlowAnswer(ans.Max, ans.Assignment);
    }

	void Stats() const {
		std::cout << "Graph with " << n << " vertices and " << m << " edges.\n";
		std::cout << "LP with " << 2 * (m + n - 2) << " constraints and " << m << " variables.\n";
	}
};

void Running(std::string test_name, std::string filename) {
	std::ifstream fin("../tests/test_data/flow/" + filename + ".ans");
    if (!fin) {
        std::cerr << "Unable to open file: " << ("../tests/test_data/flow/" + filename + ".ans") << ".\n";
        exit(EXIT_FAILURE);
    }
	long long ans;
	fin >> ans;
	fin.close();
	fin.open("../tests/test_data/flow/" + filename + ".in");
	if (!fin) {
        std::cerr << "Unable to open file: " << ("../tests/test_data/flow/" + filename + ".in") << ".\n";
        exit(EXIT_FAILURE);
    }
	int n, m, s, t;
	fin >> n >> m >> s >> t;
	std::vector<Edge> G;
	for (int i = 0; i < m; i++) {
		int from, to;
		double cap;
		fin >> from >> to >> cap;
		G.push_back(Edge(from, to, cap));
	}
	Flow flow(n, s, t, G);
	std::cout << "Flow " + test_name + " test: \n";
	flow.Stats();
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	long long answer = static_cast<long long>(flow.Solve().Max);
	const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - begin).count();
	std::cout << "Runtime (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
	std::cout << "Max: " << answer << '\n';
	if (abs(answer - ans) / ans < 1e-6) {
		std::cout << "Flow easy test passed.\n";
	} else {
		std::cout << "Flow easy test failed.\n";
	}
	fin.close();
	std::cout << "\n=============================================\n\n";
}

int main() {
	Running("easy", "flow_easy");
	Running("medium", "flow_medium");
	return 0;
}