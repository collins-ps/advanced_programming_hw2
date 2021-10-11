# homework-2
Parallelize nbody.c for a multicore processor using two different
approaches: PThreads, and OpenMP. You can determine the type of
scheduling and all other aspects of the parallelization, but the basic
idea is to have each thread (processor core) compute some fraction of
the force calculation in parallel. Getting the correct answer is a key
part of getting credit for this assignment, so be sure that you
develop a good testing strategy early. Speedup is another key aspect
-- report on performance improvements as a function of core count for
each strategy. For ease of grading, submit just two source code files
-- nbody_pt.c and nbody_omp.c.
