#include "spins.h"
using spins::O2;
using spins::O2Coupling;

#include <iostream>
#include <memory>
#include <unistd.h>

#include "testing.h"

int main(int argc, char ** argv) {
   return 0;
}

DEFINE_TEST(testSums) {
   O2 a = {1,1};
   O2 b = {0,2};
   O2 c = {1,3};

   ASSERT_EQ(spins::O2Sum(a,b), c);

   ASSERT_EQ(spins::O2Diff(c, b), a);
   ASSERT_EQ(spins::O2Diff(b, b), spins::O2_zero);
}

DEFINE_TEST(testProducts) {
   O2 a = {1,0};
   O2 b = {0,1};
   O2 c = {1,1};

   ASSERT_EQ(spins::O2Product(a,b), 0);
   ASSERT_EQ(spins::O2Product(a,c), 1);
   ASSERT_EQ(spins::O2Product(b,c), 1);
   ASSERT_EQ(spins::O2Product(c,c), 2);

   O2Coupling J = {1,2, 3,4};
   ASSERT_EQ(spins::O2Product(a,J,c), 3);
   ASSERT_EQ(spins::O2Product(b,J,c), 7);

   O2Coupling rotation = spins::O2Coupling_rotate90;
   ASSERT_EQ(spins::O2Product(a, rotation, a), 0);
   ASSERT_EQ(spins::O2Product(a, rotation, spins::O2Product(rotation, a)), -1);
   ASSERT_EQ(spins::O2_zero, spins::O2Sum(
      spins::O2Product(rotation, spins::O2Product(rotation, a)), a));
}

DEFINE_TEST(testAccumulation) {
   std::unique_ptr<O2[]> s_arr(new O2[10]);
   for (int i = 0; i < 10; ++i) {
      s_arr[i] = {(double)i, (double)2*i};
   }

   O2 sum = {0,0};
   for (int i = 0; i < 10; ++i) {
      spins::O2SumAccum(s_arr[i], &sum);
   }

   O2 answer = {45, 90};
   ASSERT_EQ(sum, answer);

   std::unique_ptr<O2Coupling[]> c_arr(new O2Coupling[10]);
   for (int i = 0; i < 5; ++i) {
      c_arr[i] = {(double)(i+1),0,0,(double)2*(i+1)};
   }
   for (int i = 0; i < 5; ++i) {
      c_arr[i+5] = {1.0/(i+1), 0, 0, 0.5/(i+1)};
   }

   O2 result = {3,4};
   answer = result;
   for (int i = 0; i < 10; ++i) {
      spins::O2ProductAccum(c_arr[i], &result);
   }

   ASSERT_EQ(result, answer);
}

DEFINE_TEST(testReflection) {
   O2 s = {21, 21};
   O2 neg_i = {-1, 0};
   ASSERT_EQ(spins::O2Reflect(spins::O2_i, spins::O2_j), neg_i);
   ASSERT_EQ(spins::O2Reflect(neg_i, neg_i), neg_i);
   ASSERT_EQ(spins::O2Reflect(spins::O2_i, s), spins::O2_j);
}
