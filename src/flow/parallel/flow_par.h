#pragma once
#include <vector>
#include <iostream>

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
    int pid, nproc;

    int n, m, s, t;
    std::vector<Edge> G;
public:
    Flow(int n, int s, int t, std::vector<Edge> G);

    FlowAnswer Solve();

	void Stats() const {
		std::cout << "Graph with " << n << " vertices and " << m << " edges.\n";
		std::cout << "LP with " << 2 * (m + n - 2) << " constraints and " << m << " variables.\n";
	}
};