#include "xychain-cli-mpi.h"

#include <gflags/gflags.h>
#include <map>
using std::map;
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <utility>
using std::pair;

#include "viewer.h"
#include "xychain-mpi.h"

bool VerifyChainArgs(const ChainArgs& args) {
   if (args.array_dims.empty() || args.island_dims.empty()) {
      return false;
   }
   if (args.n_threads < 1) {
      return false;
   }
   if (args.print_freq < 1) {
      return false;
   }
   return true;
}

bool VerifyAnnealArgs(const AnnealArgs& args) {
   if (!VerifyChainArgs(args.chain)) {
      return false;
   }

   if (args.temp_low > args.temp_high ||
       args.temp_low < 0 || args.temp_high < 0) {
      return false;
   }
   if (args.mc_sweeps < 1) {
      return false;
   }
   if (args.annealing_steps < 1) {
      return false;
   }
   if (args.overrelax_freq < 0) {
      return false;
   }
   if (args.visuals) {
      if (args.window_w < args.chain.array_dims[0]) {
         return false;
      }
      if (args.chain.array_dims.size() >= 2 && args.window_h < args.chain.array_dims[1]) {
         return false;
      }
   }

   return true;
}

bool VerifyHystArgs(const HystArgs& args) {
   if (!VerifyChainArgs(args.chain)) {
      return false;
   }

   if (args.temp <= 0) {
      return false;
   }
   if (args.field_amp == 0) {
      return false;
   }
   if (args.field_steps < 1 || args.field_cycles < 1) {
      return false;
   }
   if (args.mc_sweeps < 1) {
      return false;
   }
   if (args.visuals) {
      if (args.window_w < args.chain.array_dims[0]) {
         return false;
      }
      if (args.chain.array_dims.size() >= 2 && args.window_h < args.chain.array_dims[1]) {
         return false;
      }
   }

   return true;
}

bool VerifyPTArgs(const PTArgs& args) {
   if (!VerifyChainArgs(args.chain)) {
      return false;
   }

   if (args.temp_low > args.temp_high ||
       args.temp_low < 0 || args.temp_high < 0) {
      return false;
   }
   if (args.n_pt_replicas < 1) {
      return false;
   }
   if (args.pt_ind < 0 || args.pt_ind >= args.n_pt_replicas) {
      return false;
   }

   return true;
}

struct PTSwapComm {
   bool set_beta;
   double beta;

   int left_ind;
   int right_ind;
};

int RunPTempering(const PTArgs& pt_args) {
   if (!VerifyPTArgs(pt_args)) {
      return 1;
   }

   const ChainArgs& args = pt_args.chain;

   // Map pt-slice indices to mpi comm channels
   std::map<int, int> comm_map;
   XYChainMPI mychain(pt_args.chain.array_dims, pt_args.chain.island_dims,
                      args.array_type, args.island_type);
   
   // TODO
   return 1;
}

int RunAnnealing(const AnnealArgs& a_args) {
   if (!VerifyAnnealArgs(a_args)) {
      return 1;
   }

   const ChainArgs& args = a_args.chain;

   const int id = MPI::COMM_WORLD.Get_rank();
   const int n_proc = MPI::COMM_WORLD.Get_size();

   XYChainMPI mychain(args.array_dims, args.island_dims, args.array_type, args.island_type);
   double step_size = (a_args.temp_high - a_args.temp_low) / a_args.annealing_steps;

   if (id != 0) {
      // Slave

      mychain.Init(args.initial_config);
   
      double coupl_j = args.inter_island_coupling_amp;
      double coupl_J = args.intra_island_coupling_amp;
      if (args.has_inter_island_disorder_param) {
         mychain.InitArrayCoupling(coupl_j, args.inter_island_disorder,
                                   args.inter_island_disorder_param);
      } else {
         mychain.InitArrayCoupling(coupl_j, args.inter_island_disorder);
      }
      if (args.has_intra_island_disorder_param) {
         mychain.InitIslandCoupling(coupl_J, args.intra_island_disorder,
                                    args.intra_island_disorder_param);
      } else {
         mychain.InitIslandCoupling(coupl_J, args.intra_island_disorder);
      }
   
      mychain.InitArrayField(a_args.arr_field, a_args.arr_field_type);
      mychain.InitIslandField(a_args.isl_field, a_args.isl_field_type);

      mychain.ResyncMPI();
   
      int V = mychain.GetTotalVolume();
   
      double curtemp = a_args.temp_high;
      for (int i = 0; i < a_args.annealing_steps; ++i) {
         mychain.SetBeta(1.0/curtemp);
         if (args.proximity_form) {
            double coupl_j_new = args.inter_island_coupling_amp * 
                  exp(-args.island_spacing * sqrt(curtemp));
            double coupl_J_new = args.intra_island_coupling_amp * exp(-sqrt(curtemp));
            if (coupl_j_new != 0 && coupl_j != 0)
               mychain.ScaleArrayCouplingAmplitude(coupl_j_new / coupl_j);
            if (coupl_J_new != 0 && coupl_J != 0)
               mychain.ScaleIslandCouplingAmplitude(coupl_J_new / coupl_J);
            coupl_j = coupl_j_new;
            coupl_J = coupl_J_new;
         }
         for (int sweep = 0; sweep < a_args.mc_sweeps; ++sweep) {
            if (sweep % a_args.overrelax_freq == 0)
               mychain.MetroSweepMPI();
            else
               mychain.OverRelaxSweepMPI();
         }
         curtemp -= step_size;
      }
   } else {
      // Master
   
      XYChainMPI::Measures meas = XYChainMPI::MasterCollect(1);

      int V = mychain.GetTotalVolume();
   
      double curtemp = a_args.temp_high;
      for (int i = 0; i < a_args.annealing_steps; ++i) {
         for (int sweep = 0; sweep < a_args.mc_sweeps; ++sweep) {
            meas.m_isl_avg = 0;
            meas.m2_isl_avg = 0;
            XYChainMPI::AddMeasures(XYChainMPI::MasterCollect(1), &meas);

            double e = meas.energy / V;
            spins::O2 m = spins::O2Scale(1.0/V, meas.total_spin);
            double m2_norm = spins::O2Product(m, m);
            double m_norm = sqrt(m2_norm);
   
            double m_norm_isl_avg = meas.m_isl_avg / mychain.GetIslandVolume();
            double m2_norm_isl_avg = meas.m2_isl_avg / mychain.GetIslandVolume()
                                                               / mychain.GetIslandVolume();
   
            if (sweep % args.print_freq == 0) {
               if (args.output_type == ChainArgs::OUTPUT_BINARY) {
                  fwrite(&i, sizeof(int), 1, stdout);
                  fwrite(&sweep, sizeof(int), 1, stdout);
                  fwrite(&curtemp, sizeof(double), 1, stdout);
                  fwrite(&e, sizeof(double), 1, stdout);
                  fwrite(&m_norm, sizeof(double), 1, stdout);
                  fwrite(&m2_norm, sizeof(double), 1, stdout);
                  fwrite(&m_norm_isl_avg, sizeof(double), 1, stdout);
                  fwrite(&m2_norm_isl_avg, sizeof(double), 1, stdout);
                  fwrite(&m.x, sizeof(double), 1, stdout);
                  fwrite(&m.y, sizeof(double), 1, stdout);
               } else {
                  printf("%i %i %f ", i, sweep, curtemp);
                  printf("%f %f %f ", e, m_norm, m2_norm);
                  printf("%f %f ", m_norm_isl_avg, m2_norm_isl_avg);
                  printf("%f %f ", m.x, m.y);
                  printf("\n", m_norm_isl_avg, m2_norm_isl_avg);
               }
            }
         }
         curtemp -= step_size;
      }
   }
   MPI::Finalize();
   return 0;
}

int RunHysteresis(const HystArgs& h_args) {
   if (!VerifyHystArgs(h_args)) {
      return 1;
   }

   const ChainArgs& args = h_args.chain;

   XYChainMPI mychain(args.array_dims, args.island_dims, args.array_type, args.island_type);
   mychain.Init(args.initial_config);

   double coupl_j = args.inter_island_coupling_amp;
   double coupl_J = args.intra_island_coupling_amp;
   if (args.proximity_form) {
      coupl_j = args.inter_island_coupling_amp * 
            exp(-args.island_spacing * sqrt(h_args.temp));
      coupl_J = args.intra_island_coupling_amp * exp(-sqrt(h_args.temp));
   }

   if (args.has_inter_island_disorder_param) {
      mychain.InitArrayCoupling(coupl_j, args.inter_island_disorder,
                                args.inter_island_disorder_param);
   } else {
      mychain.InitArrayCoupling(coupl_j, args.inter_island_disorder);
   }
   if (args.has_intra_island_disorder_param) {
      mychain.InitIslandCoupling(coupl_J, args.intra_island_disorder,
                                 args.intra_island_disorder_param);
   } else {
      mychain.InitIslandCoupling(coupl_J, args.intra_island_disorder);
   }

   mychain.SetBeta(1.0/h_args.temp);

   int V = mychain.GetTotalVolume();

   double pi = acos(0)*2;
   for (int cycle = 0; cycle < h_args.field_cycles; ++cycle) {
      for (int step = 0; step < h_args.field_steps; ++step) {
         double curfield = h_args.field_amp * cos(2 * pi * step / h_args.field_steps);
         mychain.InitArrayField(curfield, XYChainMPI::CONFIG_UNIFORM);
         for (int sweep = 0; sweep < h_args.mc_sweeps; ++sweep) {
            mychain.MetroSweep();

            XYChainMPI::Measures meas = mychain.GetMeasuresExhaustive();
            double e = meas.energy / V;
            spins::O2 m = spins::O2Scale(1.0/V, meas.total_spin);
            double m2_norm = spins::O2Product(m, m);
            double m_norm = sqrt(m2_norm);
   
            double m_norm_isl_avg = meas.m_isl_avg / mychain.GetIslandVolume();
            double m2_norm_isl_avg = meas.m2_isl_avg / mychain.GetIslandVolume()
                                                                     / mychain.GetIslandVolume();
   
            if (sweep % args.print_freq == 0) {
               if (args.output_type == ChainArgs::OUTPUT_BINARY) {
                  fwrite(&cycle, sizeof(int), 1, stdout);
                  fwrite(&step, sizeof(int), 1, stdout);
                  fwrite(&sweep, sizeof(int), 1, stdout);
                  fwrite(&curfield, sizeof(double), 1, stdout);
                  fwrite(&e, sizeof(double), 1, stdout);
                  fwrite(&m_norm, sizeof(double), 1, stdout);
                  fwrite(&m2_norm, sizeof(double), 1, stdout);
                  fwrite(&m_norm_isl_avg, sizeof(double), 1, stdout);
                  fwrite(&m2_norm_isl_avg, sizeof(double), 1, stdout);
                  fwrite(&m.x, sizeof(double), 1, stdout);
                  fwrite(&m.y, sizeof(double), 1, stdout);
               } else { // if (args.output_type == ChainArgs::OUTPUT_PLAIN)
                  printf("%i %i %i %f ", cycle, step, sweep, curfield);
                  printf("%f %f %f ", e, m_norm, m2_norm);
                  printf("%f %f ", m_norm_isl_avg, m2_norm_isl_avg);
                  printf("%f %f ", m.x, m.y);
                  printf("\n");
               }
            }
         }
      }
   }

   return 0;
}
