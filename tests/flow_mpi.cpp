#include <vector>
#include <iostream>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <math.h>
#include "flow_par.h"
#include <mpi.h>

int pid, nproc;
const double epi = 1e-7;

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
		if (pid == nproc - 1) G.push_back(Edge(from, to, cap));
	}
	Flow flow(n, s, t, G);
    if (pid == nproc - 1) {
        std::cout << "Flow " + test_name + " test: \n";
        flow.Stats();
    }
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	long long answer = static_cast<long long>(flow.Solve().Max);
	const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - begin).count();
    if (pid == nproc - 1) {
        std::cout << "Runtime (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
        std::cout << "Max: " << answer << '\n';
        if ((fabs(static_cast<double>(answer) - static_cast<double>(ans)) / static_cast<double>(ans)) < epi) {
            std::cout << "Flow " + test_name + " test passed.\n";
        } else {
            std::cout << "Flow " + test_name + " test failed.\n";
        }
        std::cout << "\n=============================================\n\n";
    }
	fin.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

	Running("easy", "flow_easy");
	Running("medium", "flow_medium");

    MPI_Finalize();
	return 0;
}