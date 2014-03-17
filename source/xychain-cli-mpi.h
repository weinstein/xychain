#ifndef XYCHAIN_CLI_H
#define XYCHAIN_CLI_H

#include "xychain-mpi.h"
#include <map>
#include <string>
#include <vector>

struct ChainArgs {
   std::vector<int> array_dims;
   XYChainMPI::LatticeType array_type;
   std::vector<int> island_dims;
   XYChainMPI::LatticeType island_type;

   int n_threads;

   XYChainMPI::Config initial_config;

   double inter_island_coupling_amp;
   double intra_island_coupling_amp;

   XYChainMPI::Disorder inter_island_disorder;
   double inter_island_disorder_param;
   bool has_inter_island_disorder_param;

   XYChainMPI::Disorder intra_island_disorder;
   double intra_island_disorder_param;
   bool has_intra_island_disorder_param;

   bool proximity_form;
   double island_spacing;

   int print_freq;

   enum OutputType {
      OUTPUT_PLAIN,
      OUTPUT_BINARY
   };
   OutputType output_type;
};

struct AnnealArgs {
   ChainArgs chain;

   double temp_low;
   double temp_high;

   double arr_field;
   XYChainMPI::Config arr_field_type;
   double isl_field;
   XYChainMPI::Config isl_field_type;

   int mc_sweeps;
   int annealing_steps;
   int overrelax_freq;

   bool visuals;
   int window_w;
   int window_h;
};

struct HystArgs {
   ChainArgs chain;

   double temp;

   double field_amp;
   int field_steps;
   int field_cycles;

   int mc_sweeps;

   bool visuals;
   int window_w;
   int window_h;
};

struct PTArgs {
   ChainArgs chain;

   double temp_low;
   double temp_high;

   int mc_sweeps;
   int overrelax_freq;

   int n_pt_replicas;
   int pt_ind;
};

class Command {
 public:
   Command() { }
   virtual int Run() = 0;
   virtual int RunHelp() = 0;
};

extern std::map<std::string, Command*> xychain_commands;

XYChainMPI::Config ParseXYChainConfig(const std::string& str);
XYChainMPI::Disorder ParseXYChainDisorder(const std::string& str);

bool VerifyChainArgs(const ChainArgs& args);
bool VerifyAnnealArgs(const AnnealArgs& args);
bool VerifyHystArgs(const HystArgs& args);
bool VerifyPTArgs(const PTArgs& args);
int RunAnnealing(const AnnealArgs& args);
int RunHysteresis(const HystArgs& args);
int RunPTempering(const PTArgs& args);

#endif
