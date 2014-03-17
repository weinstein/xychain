#include "rng.h"
using rng::RandomGen;

#include <iostream>
#include <unistd.h>
using namespace std;
#include "testing.h"

int main() {
   return 0;
}

DEFINE_TEST(testManyRandoms) {
   RandomGen rng;
   int buckets[100];
   for (int i = 0; i < 100; ++i) {
      buckets[i] = 0;
   }
   for (int i = 0; i < 1e6; ++i) {
      double rand = rng.RandomUniform();
      ASSERT_GTEQ(rand, 0);
      ASSERT_LTEQ(rand, 1);
      int b_ind = (int)(rand*100.0);
      ++buckets[b_ind];
   }
   for (int i = 0; i < 100; ++i) {
      int diff = buckets[i] - 1e4;
      ASSERT_LTEQ(diff, 500);
      ASSERT_GTEQ(diff, -500);
   }
}

DEFINE_TEST(testActuallyRandom) {
   RandomGen rng;
   sleep(0.1);
   RandomGen rng2;
   ASSERT_NEQ(rng.RandomUniform(), rng2.RandomUniform());
}
