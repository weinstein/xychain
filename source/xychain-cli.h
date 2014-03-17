#ifndef XYCHAIN_CLI_H
#define XYCHAIN_CLI_H

#include "xychain.h"
#include <map>
#include <string>
#include <vector>

struct ChainArgs {
   std::vector<int> array_dims;
   std::vector<int> island_dims;

   int n_threads;

   XYChain::Config initial_config;

   double inter_island_coupling_amp;
   double intra_island_coupling_amp;

   XYChain::Disorder inter_island_disorder;
   double inter_island_disorder_param;
   bool has_inter_island_disorder_param;

   XYChain::Disorder intra_island_disorder;
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
   XYChain::Config arr_field_type;
   double isl_field;
   XYChain::Config isl_field_type;

   int mc_sweeps;
   int annealing_steps;
   int overrelax_freq;

   bool visuals;
   int window_w;
   int window_h;

   std::string png_pref;
   int frame_freq;
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

   std::string png_pref;
   int frame_freq;
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

XYChain::Config ParseXYChainConfig(const std::string& str);
XYChain::Disorder ParseXYChainDisorder(const std::string& str);

bool VerifyChainArgs(const ChainArgs& args);
bool VerifyAnnealArgs(const AnnealArgs& args);
bool VerifyHystArgs(const HystArgs& args);
bool VerifyPTArgs(const PTArgs& args);
int RunAnnealing(const AnnealArgs& args);
int RunHysteresis(const HystArgs& args);
int RunPTempering(const PTArgs& args);

#endif
