#include <vector>
#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "lp_seq.h"
#include "flow_seq.h"

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
	FlowSeq flow(n, s, t, G);
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