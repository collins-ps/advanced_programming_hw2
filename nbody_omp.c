#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include "timer.h"

#define SOFTENING 1e-9f
#define nt 8
#define MARGIN 0.0001 // margin of error allowed between parallel-calculated and serial-calculated values

typedef struct {
  float x, y, z;        /* particle positions */
  float vx, vy, vz;     /* particle momenta */
} Particle;


/* randomly initialize particle positions and momenta */
void ran_init(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

/* calculate all interparticle forces and update instantaneous velocities */
void calc_force(Particle *p, float dt, int n) {
    int i;
    #pragma omp parallel for default(shared) private(i) num_threads(nt)
    for (int i = 0; i < n; i++) { 
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      /* calculate net particle for on i'th particle */
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }
    /* update instantaneous velocity based on force and timestep */
    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

/* calculate all interparticle forces and update instantaneous velocities */
void calc_force_serial(Particle *p, float dt, int n) {
    for (int i = 0; i < n; i++) { 
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      /* calculate net particle for on i'th particle */
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }
    /* update instantaneous velocity based on force and timestep */
    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
}

int main(const int argc, const char** argv) {
  FILE *datafile    = NULL;      /* output file for particle positions */
  //FILE *resultsfile    = NULL;  /* results file for performance tracking */
  int   nParticles  = 3000;      /* number of particles */

  if (argc > 1)
    nParticles      = atoi(argv[1]);

  const float dt    = 0.01f; /* time step   */
  const int nIters  = 200;   /* number of steps in simulation */

  float *buf        =  malloc(nParticles*sizeof(Particle));
  Particle  *p          = (Particle *) buf;

  ran_init(buf, 6*nParticles); /* Init pos and vel data */

  // create copy of buf for testing purposes //
  float *buf_test        =  malloc(nParticles*sizeof(Particle));
  Particle  *p_test          = (Particle *) buf_test;
  for (int k = 0; k < 6*nParticles; k++){
      buf_test[k] = buf[k];
  }

  double totalTime  = 0.0;
  double baseTime_total = 0.0;

  datafile          = fopen("particles.dat","w");
  fprintf(datafile,"%d %d %d\n", nParticles, nIters, 0);

  //resultsfile          = fopen("results_omp.dat","a");
  // fprintf(resultsfile,"%15s %15s %15s %15s\n", "nThreads", "avgTime", "totalTime", "Performance Improvement, Percent");

  /* ------------------------------*/
  /*     MAIN LOOP                 */
  /* ------------------------------*/
  for (int iter = 1; iter <= nIters; iter++) {
    printf("iteration:%d\n", iter);
    
    for (int i = 0;i < nParticles; ++i)
      fprintf(datafile, "%f %f %f \n", p[i].x, p[i].y, p[i].z);
    
    double omp_get_wtime(void);
    double start; 
    double end; 
    start = omp_get_wtime(); 

    calc_force(p, dt, nParticles);           /* compute interparticle forces and update vel */

    for (int i = 0 ; i < nParticles; i++) {  /* compute new position */
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }
    
    end = omp_get_wtime(); 
    if (iter > 1) {                          /* First iter is warm up */
      totalTime += (end - start); 
    }

    // test performance and accuracy //
    StartTimer();
    calc_force_serial(p_test, dt, nParticles);
    for (int i = 0 ; i < nParticles; i++) {  /* compute new position */
      p_test[i].x += p_test[i].vx*dt;
      p_test[i].y += p_test[i].vy*dt;
      p_test[i].z += p_test[i].vz*dt;
    }
    const double tElapsed_base = GetTimer() / 1000.0;
      if (iter > 1) {                          
        baseTime_total += tElapsed_base; 
      } 
    for (int i = 0; i < nParticles; i++){
        assert(fabs(p[i].vx - p_test[i].vx) < MARGIN); 
        assert(fabs(p[i].vy - p_test[i].vy) < MARGIN); 
        assert(fabs(p[i].vz - p_test[i].vz) < MARGIN); 
    }
    printf("Passed accuracy test for iteration %d.\n",iter); 

  }

  // test accuracy //
  for (int i = 0 ; i < nParticles; i++) { 
      assert(fabs(p[i].vx - p_test[i].vx) < MARGIN); 
      assert(fabs(p[i].vy - p_test[i].vy) < MARGIN); 
      assert(fabs(p[i].vz - p_test[i].vz) < MARGIN); 
    }
  printf("Passed all tests.\n"); 
  
  fclose(datafile);
  double avgTime = totalTime / (double)(nIters-1); 
  double avgTime_base = baseTime_total / (double)(nIters-1);
  double perf_imprv = fabs(avgTime_base - avgTime)/ avgTime * 100;

  printf("avgTime: %f   totTime: %f \n", avgTime, totalTime);
  //fprintf(resultsfile, "%15d %15f %15f %15.2f \n", nt, avgTime, totalTime, perf_imprv);
  //fclose(resultsfile);
  free(buf);
  free(buf_test);
}