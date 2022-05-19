/******************************************************************************  
 *                                                                            *  
 * RANDOM.C                                                                   *  
 *                                                                            *  
 * WRAPPERS FOR RANDOM NUMBER GENERATOR                                       *  
 *                                                                            *  
 ******************************************************************************/

#include "decs.h"

static gsl_rng **rng;

// Use Mersenne twister

/**
 * @brief Initializes a differently seeded Mersenne Twister generator on every thread
 * 
 * @param seed each thread is seeded with (seed + thread number)
 */
void init_random(int seed)
{
  rng = malloc(nthreads*sizeof(gsl_rng*));
  #pragma omp parallel
  {
    rng[omp_get_thread_num()] = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng[omp_get_thread_num()], seed + omp_get_thread_num());
  }
}

double get_rand()
{
  return gsl_rng_uniform(rng[omp_get_thread_num()]);
}

double get_chisq(double nu)
{
  return gsl_ran_chisq(rng[omp_get_thread_num()], nu);
}

void get_ran_dir_3d(double *nx, double *ny, double *nz)
{
  gsl_ran_dir_3d(rng[omp_get_thread_num()], nx, ny, nz);
}

