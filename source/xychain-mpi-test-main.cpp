#include <iostream>
using std::cout;
#include <mpi.h>
#include <string>
using std::string;
#include <vector>
using std::vector;

#include "testing.h"
#include "xychain-mpi.h"

void RunTest();

int main(int argc, char** argv) {
   MPI::Init(argc, argv);

   //RunTest();

   const int id = MPI::COMM_WORLD.Get_rank();
   const int n_proc = MPI::COMM_WORLD.Get_size();

   vector<int> arr_dims = {16,16};
   vector<int> isl_dims = {1};
   if (id == 0) {
      // I am the master
      XYChainMPI::Measures measures;
      XYChainMPI::ZeroMeasures(&measures);

      cout << "Resync\n";
      XYChainMPI::AddMeasures(XYChainMPI::MasterResync(), &measures);
      cout << "Sweep(3000)\n";
      XYChainMPI::AddMeasures(XYChainMPI::MasterMetroSweeps(3000), &measures);

      cout << "Sweeping..\n";
      double m2_avg = 0;
      double e_avg = 0;
      const int V = arr_dims[0]*arr_dims[1];
      const int n_sweep = 1000;
      for (int i = 0; i < n_sweep; ++i) {
         measures.m_isl_avg = 0;
         measures.m2_isl_avg = 0;
         XYChainMPI::AddMeasures(XYChainMPI::MasterMetroSweeps(1), &measures);
         e_avg += measures.energy / V;
         spins::O2 m = spins::O2Scale(1.0 / V, measures.total_spin);
         m2_avg += spins::O2Product(m, m);
         ASSERT_LT(measures.m2_isl_avg, 1.000001);
         ASSERT_GT(measures.m2_isl_avg, 0.999999);
      }
      e_avg = e_avg / n_sweep;
      m2_avg = m2_avg / n_sweep;
      ASSERT_GT(e_avg, -0.6);
      ASSERT_LT(e_avg, -0.5);
      ASSERT_LT(m2_avg, 0.05);

      cout << "Shutdown\n";
      XYChainMPI::MasterShutdownSlaves();
      cout << "Success!\n";
   } else {
      const int my_id = MPI::COMM_WORLD.Get_rank();

      // I am a slave
      XYChainMPI mychain(arr_dims, isl_dims);

      mychain.Init(XYChainMPI::CONFIG_RANDOM);
      mychain.ScaleArrayCouplingAmplitude(1);
      mychain.ScaleIslandCouplingAmplitude(0);
      mychain.SetBeta(0.5);  // below critical beta

      mychain.CommandLoop();
      fprintf(stderr, "(%i): finished\n", my_id);
   }
}

void RunTest() {
   vector<int> dims16x16 = {16, 16};
   vector<int> dims8x10x12 = {8, 10, 12};
   XYChainMPI chain_a(dims16x16, dims8x10x12);
   XYChainMPI chain_b(dims16x16, dims8x10x12);

   chain_a.Init(XYChainMPI::CONFIG_RANDOM);
   chain_a.InitArrayCoupling(1, XYChainMPI::DISORDER_UNIFORM_PHASE);
   chain_a.InitIslandCoupling(1, XYChainMPI::DISORDER_UNIFORM_PHASE);
   
   chain_b.Init(XYChainMPI::CONFIG_RANDOM);
   chain_b.InitArrayCoupling(1, XYChainMPI::DISORDER_UNIFORM_PHASE);
   chain_b.InitIslandCoupling(1, XYChainMPI::DISORDER_UNIFORM_PHASE);

   const int arr_coord = 4;
   const int isl_coord = 6;
   const int n_arr_slice = chain_a.GetNarr();
   const int n_spins = n_arr_slice * chain_a.GetIslandVolume();
   const int buflen = 1 + 2 * sizeof(int)
                        + n_arr_slice * sizeof(spins::O2)
                        + n_arr_slice * arr_coord * sizeof(spins::O2Coupling)
                        + n_arr_slice * arr_coord * sizeof(int)
                        + n_spins * sizeof(spins::O2)
                        + n_spins * isl_coord * sizeof(spins::O2Coupling)
                        + n_spins * isl_coord * sizeof(int);
   //cout << "BTW, buflen=" << buflen << endl;
   char* buf = new char[buflen];
   buf[buflen-1] = 0;
   chain_a.FlattenStaticInfo(0, n_arr_slice/2, buf);
   chain_b.UnflattenStaticInfo(buf);
   chain_b.FlattenStaticInfo(n_arr_slice/2, n_arr_slice, buf);
   chain_a.UnflattenStaticInfo(buf);
   delete [] buf;

   for (int i = 0; i < chain_a.GetNarr(); ++i) {
      // Array neighboring islands
      for (int j = 0; j < arr_coord; ++j) {
         ASSERT_EQ(chain_a.GetArrayNeighborIndex(i,j),
                     chain_b.GetArrayNeighborIndex(i,j));
         ASSERT_EQ(chain_a.GetArrayNeighborCoupling(i,j),
                     chain_b.GetArrayNeighborCoupling(i,j));
         
         int neighbor_ind = chain_a.GetArrayNeighborIndex(i, j);
         int alt_j = (j + arr_coord/2) % (arr_coord);
         ASSERT_EQ(chain_a.GetArrayNeighborCoupling(i, j),
                   spins::O2CouplingTranspose(
               chain_a.GetArrayNeighborCoupling(neighbor_ind, alt_j)));
         ASSERT_EQ(chain_b.GetArrayNeighborCoupling(i, j),
                   spins::O2CouplingTranspose(
               chain_b.GetArrayNeighborCoupling(neighbor_ind, alt_j)));
      }

      for (int site = 0; site < chain_a.GetIslandVolume(); ++site) {
         // Island neighboring sites
         for (int j = 0; j < isl_coord; ++j) {
            ASSERT_EQ(chain_a.GetIslandSiteNeighborIndex(i, site, j),
                        chain_b.GetIslandSiteNeighborIndex(i, site, j));
            ASSERT_EQ(chain_a.GetIslandSiteNeighborCoupling(i, site, j),
                        chain_b.GetIslandSiteNeighborCoupling(i, site, j));

            int neighbor_site = chain_a.GetIslandSiteNeighborIndex(i, site, j);
            int alt_j = (j + isl_coord/2) % (isl_coord);
            ASSERT_EQ(chain_a.GetIslandSiteNeighborCoupling(i, site, j),
                      spins::O2CouplingTranspose(
                  chain_a.GetIslandSiteNeighborCoupling(i, neighbor_site, alt_j)));
            ASSERT_EQ(chain_b.GetIslandSiteNeighborCoupling(i, site, j),
                      spins::O2CouplingTranspose(
                  chain_b.GetIslandSiteNeighborCoupling(i, neighbor_site, alt_j)));
         }
      }
   }
}
