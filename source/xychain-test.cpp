#include "xychain.h"

#include <iostream>
#include "testing.h"
#include <vector>
using std::vector;

int main() {
   return 0;
}

DEFINE_TEST(testChainInit) {
   vector<int> arr_dims = {32, 32};
   vector<int> isl_dims = {1};
   XYChain mychain(arr_dims, isl_dims, 1);

   ASSERT_EQ(mychain.GetTotalVolume(), 32*32);
   ASSERT_EQ(mychain.GetIslandVolume(), 1);
   ASSERT_EQ(mychain.GetNarr(), 32*32);

   mychain.Init(XYChain::CONFIG_UNIFORM);
   mychain.ScaleArrayCouplingAmplitude(1);
   mychain.ScaleIslandCouplingAmplitude(0);

   XYChain::Measures isl0_meas = mychain.GetIslandMeasures(0);
   ASSERT_EQ(isl0_meas.total_spin, spins::O2_i);
   ASSERT_EQ(mychain.GetMeasures().total_spin,
             spins::O2Scale(32*32, spins::O2_i));
   ASSERT_EQ(mychain.GetMeasures().energy,
             -2*32*32);

   vector<int> arr_dims_2 = {1, 1, 1};
   vector<int> isl_dims_2 = {16, 64};
   XYChain mychain2(arr_dims_2, isl_dims_2, 1);

   ASSERT_EQ(mychain2.GetTotalVolume(), 16*64);
   ASSERT_EQ(mychain2.GetIslandVolume(), 16*64);
   ASSERT_EQ(mychain2.GetNarr(), 1);

   mychain2.Init(XYChain::CONFIG_UNIFORM);
   mychain2.ScaleArrayCouplingAmplitude(0);
   mychain2.ScaleIslandCouplingAmplitude(1);

   XYChain::Measures isl_meas_2 = mychain2.GetIslandMeasures(0);
   ASSERT_EQ(isl_meas_2.total_spin, spins::O2Scale(16*64, spins::O2_i));
   ASSERT_EQ(mychain2.GetMeasures().total_spin,
             spins::O2Scale(16*64, spins::O2_i));
   ASSERT_EQ(mychain.GetMeasures().energy,
             -2*16*64);

   XYChain mychain3(arr_dims_2, isl_dims_2, 1);
   mychain3.Init(XYChain::CONFIG_RANDOM);
   spins::O2 m = mychain3.GetMeasuresExhaustive().total_spin;
   spins::O2ScaleAccum(1.0 / mychain3.GetTotalVolume(), &m);
   double total_spin_mag = spins::O2Product(m, m);
   ASSERT_LT(total_spin_mag, 0.01);
}

DEFINE_TEST(testNeighborInitializationSquare) {
   vector<int> arr_dims = {3, 3};
   vector<int> isl_dims = {1};
   XYChain mychain(arr_dims, isl_dims, 1);

   ASSERT_EQ(mychain.GetArrayNeighborIndex(1, 3), 7);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 1), 3);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(8, 0), 6);
}

DEFINE_TEST(testNeighborInitializationRect) {
   vector<int> arr_dims = {2, 4, 3};
   vector<int> isl_dims = {1};
   XYChain mychain(arr_dims, isl_dims, 1);

   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 0), 1);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 1), 2);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 2), 8);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 3), 1);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 4), 6);
   ASSERT_EQ(mychain.GetArrayNeighborIndex(0, 5), 16);
}

DEFINE_TEST(testMCSweep) {
   vector<int> arr_dims = {16, 16};
   vector<int> isl_dims = {1};
   XYChain mychain(arr_dims, isl_dims, 4);

   mychain.Init(XYChain::CONFIG_RANDOM);
   mychain.ScaleArrayCouplingAmplitude(1);
   mychain.ScaleIslandCouplingAmplitude(0);
   mychain.SetBeta(0.5);  // below critical beta

   ASSERT_EQ(mychain.GetTotalVolume(), 16*16);

   const int eq_sweep = 2000;
   const int n_sweep = 1000;
   int V = mychain.GetTotalVolume();
   double m2_avg = 0;
   double e_avg = 0;
   for (int i = 0; i < eq_sweep; ++i) mychain.MetroSweep();
   for (int i = 0; i < n_sweep; ++i) {
      mychain.MetroSweep();
      e_avg += mychain.GetMeasures().energy / V;
      spins::O2 m = spins::O2Scale(1.0 / V, mychain.GetMeasures().total_spin);
      m2_avg += spins::O2Product(m, m);
   }
   ASSERT_EQ((int)mychain.GetMeasures().energy,
             (int)mychain.GetMeasuresExhaustive().energy);
   e_avg = e_avg / n_sweep;
   m2_avg = m2_avg / n_sweep;
   ASSERT_GT(e_avg, -0.6);
   ASSERT_LT(e_avg, -0.5); 
   ASSERT_LT(m2_avg, 0.05);

   mychain.SetBeta(2.0);  // nearer to critical beta
   e_avg = 0;
   for (int i = 0; i < eq_sweep; ++i) mychain.MetroSweep();
   for (int i = 0; i < n_sweep; ++i) {
      mychain.MetroSweep();
      e_avg += mychain.GetMeasures().energy / V;
   }
   ASSERT_EQ((int)mychain.GetMeasures().energy,
             (int)mychain.GetMeasuresExhaustive().energy);
   e_avg = e_avg / n_sweep;
   ASSERT_GT(e_avg, -1.8);
   ASSERT_LT(e_avg, -1.6);
}

DEFINE_TEST(testDisorderInit) {
   vector<int> dims16x16 = {16, 16};
   vector<int> dims8x10x12 = {8, 10, 12};
   XYChain mychain(dims16x16, dims8x10x12, 1);

   mychain.Init(XYChain::CONFIG_RANDOM);
   mychain.InitArrayCoupling(1, XYChain::DISORDER_UNIFORM_PHASE);
   mychain.InitIslandCoupling(1, XYChain::DISORDER_UNIFORM_PHASE);

   int arr_coord = 4;
   int isl_coord = 6;
   for (int i = 0; i < mychain.GetNarr(); ++i) {
      // Array neighboring islands
      for (int j = 0; j < arr_coord; ++j) {
         int neighbor_ind = mychain.GetArrayNeighborIndex(i, j);
         int alt_j = (j + arr_coord/2) % (arr_coord);
         ASSERT_EQ(mychain.GetArrayNeighborCoupling(i, j),
                   spins::O2CouplingTranspose(
               mychain.GetArrayNeighborCoupling(neighbor_ind, alt_j)));
         spins::O2Coupling c = mychain.GetArrayNeighborCoupling(i, j);
         ASSERT_EQ((int)(1e6 * (c.a * c.d - c.b * c.c) + 0.5), 1e6);
      }

      for (int site = 0; site < mychain.GetIslandVolume(); ++site) {
         // Island neighboring sites
         for (int j = 0; j < isl_coord; ++j) {
            int neighbor_site = mychain.GetIslandSiteNeighborIndex(i, site, j);
            int alt_j = (j + isl_coord/2) % (isl_coord);
            ASSERT_EQ(mychain.GetIslandSiteNeighborCoupling(i, site, j),
                      spins::O2CouplingTranspose(
                  mychain.GetIslandSiteNeighborCoupling(i, neighbor_site, alt_j)));
            spins::O2Coupling c = mychain.GetIslandSiteNeighborCoupling(i, site, j);
            ASSERT_EQ((int)(1e6 * (c.a * c.d - c.b * c.c) + 0.5), 1e6);
         }
      }
   }

   ASSERT_EQ((int)mychain.GetMeasures().energy,
             (int)mychain.GetMeasuresExhaustive().energy);
   for (int i = 0; i < 10; ++i) {
      mychain.MetroSweep();
      ASSERT_EQ((int)mychain.GetMeasures().energy,
                (int)mychain.GetMeasuresExhaustive().energy);
   }
   for (int i = 0; i < 10; ++i) {
      mychain.OverRelaxSweep();
      ASSERT_EQ((int)mychain.GetMeasures().energy,
                (int)mychain.GetMeasuresExhaustive().energy);
   }
}
