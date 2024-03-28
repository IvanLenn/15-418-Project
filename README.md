<p align="center">

  <h1 align="center">Parallel Linear Programming Solver Using MPI</h1>
  <p align="center">
    <strong>Yuan (Alex) Guan</strong></a>
    ·
    <strong>Yixuan (Ivan) Zheng</strong></a>
    <br>
    <a href="https://github.com/IvanLenn/15-418-Project/blob/main/report/418proposal.pdf"><strong>Project proposal</strong></a>
  <div align="center"></div>
</p>

## SUMMARY
We are going to implement a parallel version of a parallel linear programming solver with MPI. We will parallelize Simplex algorithm in order to speed up the computation of linear programming solver. We plan to experiment on the GHC Clusters and PSC machines and perform analysis on the results. We'll probably study the parallelization of other LP algorithms such as interior point algorithms and Seidal's algorithm.

## BACKGROUND: 
Linear programming(LP) is a technique for the optimization of a linear objective function. The development of better and faster LP algorithms is an interesting ongoing area of research. 

The Simplex algorithm is a common method to solve LP. Starting with a corner of the feasible region, we repeatedly look at all neighboring corners of our current position and go to the best one. When we reach a position without better neighboring corners, it means we have found the optimal solution because the feasible region is convex.
    This algorithm can be greatly customized by large numbers of parameters and tolerances. We wish to speed up the algorithm using parallel techniques with message-passing.
     
    

## THE CHALLENGE: 
- Data Partitioning: Divide the constraints and objective function coefficients among MPI processes. Use data distribution techniques to distribute the data in a balanced way and minimize communication overhead. Locality can be a problem for LP algorithms dependent on the input constraints.

- Load Balancing: Implement load balancing work assignment methods to distribute the computational workload among MPI processes. This may involve dynamically redistributing data or adjusting the workload of each process based on its progress during the algorithm's execution. We'll experiment with different levels of dynamic assignment from static to fine-grained.

- Communication Optimization: Minimize communication overhead by aggregating data exchanges between processes using techniques such as batches and reducing the frequency of artifactual communication where possible.

- Parallelization Strategy: Explore different parallelization strategies, such as parallelizing different stages of the Simplex algorithm or parallelizing multiple independent Simplex iterations concurrently. We'll probably study the parallelization of other LP algorithms such as interior point algorithms and Seidal's algorithm.

## RESOURCES: 
- Type of computers: GHC for $\leq$ 8 cores. PSC for multi-cores experiments.
- We will probably use some existing code for the serial version implementation.
- We are going to implement the parallel version with MPI.

## GOALS AND DELIVERABLES: 
- Basic Goal:
    A working serial version of the Simplex algorithm, a benchmark script to measure the performance, as well as a working serial version of the Simplex algorithm using MPI.
- Ideal Goal:
    Optimize our parallel version by using proper parallel strategy and techniques after experiments and analysis. We wish to achieve a 10x speedup on PSC machines with high thread counts. (Past research using the shared memory approach has reached 19x speedup)
- Extra Goal:
    If we've achieved the ideal goal, we'll try to further optimize our parallel version of the Simplex algorithm if possible. We'll probably utilize efficient (perhaps parallel lock-free) data structures to help us speed up the program as well.\\

## PLATFORM CHOICE: 
- We’ll use GHC machines for 8-core CPU testing for experiments for our basic implementation. 
- We’ll also use PSC machines for experiments with multi-nodes.
- We'll C++ for this task, which is also appropriate since C++ has excellent support for MPI interface for parallel programming with MPI.

## SCHEDULE: 
| Week | Plan | Status | Comment | 
|------|------|--------|---------|
| March 25 | Finalize Project Idea | :white_check_mark: | **Project Proposal Due: Wednesday, March 27th, 11:59pm** |
| April 1 | Implement both serial and parallel version, do various experiments |        |         |
| April 8|      |        |         |
| April 15 | Compare different approaches of parallelism and do scalability analysis |        | **Milestone Report Due: Tuesday, April 16th, 11:59pm** |
| April 22 | Extra optimization |        |         |
| April 29 |      |        |         |
| May 6 |      |        | **Poster Session: Monday, May 6th, 1:00-4:00pm <br> Final Report Due: Sunday, May 5th, 11:59pm** |

