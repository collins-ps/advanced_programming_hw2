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

# MPCS 51100 HW2
**Name**: Phoebe Collins, UCID: 12277438

## References
I did not use any references other than class materials.

## Discussion on program
I did not observe any compilation or run-time errors. 

In comparing the value of vx, vy, and vz generated by the parallel calc_force function and those generated by the serial calc_force function, I allow for a margin of error of up to 0.0001. This is to account for the non-associativity of floating point arithmetic and the chaos of the n-body problem which means that even small differences in the values could grow exponentially with time. 

## Compiling and running
`make nbody_pt` and then `./nbody_pt` which uses the default 3000 particles or `./nbody_pt x` where x is an integer representing the number of particles. Finally, run `make clean`.

`make nbody_omp` and then `./nbody_omp` which uses the default 3000 particles or `./nbody_omp x` where x is an integer representing the number of particles. Finally, run `make clean`.

## Performance comparison for PThreads strategy
Program was run on MPCS Linux clusters without any optimization flags to derive these timings. Performance improvement indicates the percentage decrease in average time taken by the parallel program over 200 iterations as compared to the average time taken by the original serial program for the same number of particles, number of iterations, and same initial pos and vel data.

          nThreads     avgTime        totalTime    Performance Improvement(%)
              8        0.034928        6.950679          302.30 
              4        0.060246       11.988948          135.04 
              2        0.107256       21.343989           33.76 
              1        0.151393       30.127192            5.86 
              
## Performance comparison for OpenMP strategy
Program was run on MPCS Linux clusters without any optimization flags to derive these timings. Performance improvement indicates the percentage decrease in average time taken by the parallel program over 200 iterations as compared to the average time taken by the original serial program for the same number of particles, number of iterations, and same initial pos and vel data.

          nThreads      avgTime       totalTime    Performance Improvement(%)
              8        0.033829        6.732069          317.21 
              4        0.041190        8.196903          239.19 
              2        0.073572       14.640751           84.54 
              1        0.138812       27.623505            0.05 
