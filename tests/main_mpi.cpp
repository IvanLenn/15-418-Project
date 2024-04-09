#include <iostream>
#include <vector>
#include <unistd.h>
#include <string>
#include <fstream>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int pid, nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    std::cout << "Hello from process " << pid << "\n";

    MPI_Finalize();
    return 0;
}