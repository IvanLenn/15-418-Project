<p align="center">

  <h1 align="center">Parallel Linear Programming Solver Using MPI</h1>
  <p align="center">
    <strong>Yuan (Alex) Guan</strong></a>
    ·
    <strong>Yixuan (Ivan) Zheng</strong></a>
    <br>
    <a href="doc/Proposal.pdf"><strong>Project proposal</strong></a>
	<a href="doc/Final_report.pdf"><strong>Final report</strong></a>
  <div align="center"></div>
</p>


## SUMMARY
A parallel linear programming with modified Simplex algorithm using Message Passing Interface(MPI). Performance analysis are conducted on 8-core GHC machines and up to 128-core PSC machines. On general problems with 300 varables and 6000 constraints achieve over 70x speedup on 128 cores with column-wise data partitioning optimization. 


## BACKGROUND: 
Linear programming(LP) is a technique for the optimization of a linear objective function. The development of better and faster LP algorithms is an interesting ongoing area of research. 

The Simplex algorithm is a common method to solve LP. Starting with a corner of the feasible region, we repeatedly look at all neighboring corners of our current position and go to the best one. When we reach a position without better neighboring corners, it means we have found the optimal solution because the feasible region is convex.
    This algorithm can be greatly customized by large numbers of parameters and tolerances. We wish to speed up the algorithm using parallel techniques with message-passing.
    

## CHALLENGE: 
- Data Partitioning: Divide the constraints and objective function coefficients among MPI processes. Use data distribution techniques to distribute the data in a balanced way and minimize communication overhead. Locality can be a problem for LP algorithms dependent on the input constraints.

- Load Balancing: Implement load balancing work assignment methods to distribute the computational workload among MPI processes. This may involve dynamically redistributing data or adjusting the workload of each process based on its progress during the algorithm's execution. We'll experiment with different levels of dynamic assignment from static to fine-grained.

- Communication Optimization: Minimize communication overhead by aggregating data exchanges between processes using techniques such as batches and reducing the frequency of artifactual communication where possible.

- Parallelization Strategy: Explore different parallelization strategies, such as parallelizing different stages of the Simplex algorithm or parallelizing multiple independent Simplex iterations concurrently. We'll probably study the parallelization of other LP algorithms such as interior point algorithms and Seidal's algorithm.


## POSTER
<img src="doc/Poster.png" alt="Poster" width="1000">
