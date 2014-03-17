#ifndef RANDOM_H
#define RANDOM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>

namespace rng {

class RandomGen {
 public:
   RandomGen() {
      fRNG = gsl_rng_alloc(gsl_rng_mt19937);
      struct timeval t;
      gettimeofday(&t, NULL);
      unsigned long long int seed = t.tv_usec + 1000000LL * t.tv_sec;
      gsl_rng_set(fRNG, seed);
   }

   ~RandomGen() {
      gsl_rng_free(fRNG);
   }

   inline double RandomUniform() {
      return gsl_rng_uniform(fRNG);
   }

   inline double RandomNormal() {
      return gsl_ran_gaussian(fRNG, 1);
   }

 private:
   gsl_rng* fRNG;
};

}  // namespace rng

#endif
