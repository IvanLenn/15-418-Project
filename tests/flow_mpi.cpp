#include <vector>
#include <unistd.h>
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

void Running(std::string test_name, std::string filename, std::string out_file, int verbose) {
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
    if (pid == nproc - 1 && verbose) {
        std::cout << "Flow " + test_name + " test: \n";
        flow.Stats();
    }
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	long long answer = static_cast<long long>(flow.Solve().Max);
	const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - begin).count();
    if (pid == 0 && out_file != "") {
        std::ofstream fout("../tests/test_data/flow/" + out_file + ".txt", std::ios::app);
        fout << init_time << ", ";
        fout.close();
    }
    if (pid == nproc - 1 && verbose) {
        std::cout << "Runtime (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
        std::cout << "Max: " << answer << '\n';
        if ((fabs(static_cast<double>(answer) - static_cast<double>(ans)) / static_cast<double>(ans)) < epi) {
            std::cout << "Flow " + test_name + " test passed.\n";
        } else {
            std::cout << "Flow " + test_name + " test failed.\n";
        }
        std::cout << "\n=============================================\n\n";
    }
    else {
        if (pid == nproc - 1 && (fabs(static_cast<double>(answer) - static_cast<double>(ans)) / static_cast<double>(ans)) > epi) {
            std::cout << "Flow " + test_name + " test failed.\n";
        }
    }
	fin.close();
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    int opt, verbose = 1;
    std::string file_name = "";
    int running_tag = 7;
	while ((opt = getopt(argc, argv, "v:f:t:")) != -1) {
		switch (opt) {
		case 'v':
			verbose = atoi(optarg);
			break;
		case 'f':
            file_name = optarg;
            break;
        case 't':
            running_tag = atoi(optarg);
            break;
		default:
			if (pid == 0) {
				std::cerr << "Usage: " << argv[0] << " [-v verbose] [-f out_file] [-t running tag]\n";
			}

			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
	}

    std::vector<std::string> test_names = {"easy", "medium", "hard"};
    std::vector<std::string> file_names = {"flow_easy", "flow_medium", "flow_hard"};
    for (int i = 0; i < 3; i++) {
        if (running_tag & (1 << i)) {
            Running(test_names[i], file_names[i], file_name, verbose);
        }
    }

    MPI_Finalize();
	return 0;
}