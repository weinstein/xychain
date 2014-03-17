#include "xychain-cli-mpi.h"
#include "parsestring.h"
#include <gflags/gflags.h>

#include <stdio.h>
#include <string>
using std::string;
#include <map>
using std::map;
#include <math.h>
#include <vector>
using std::vector;
#include <utility>
using std::pair;

#include <mpi.h>

// Define ALL the flags!
// XYChain parameters, which are independent of the command being run
DEFINE_string(array_dims, "", "comma-seperated list of array dimensions");
DEFINE_string(island_dims, "", "comma-seperated list of island dimensions");
DEFINE_double(inter_island_coupling, 1, "coupling amplitude between islands");
DEFINE_double(intra_island_coupling, 1, "coupling amplitude between spins on an island");
DEFINE_string(initial_config, "random", "initial configuration type: uniform, random, or a file "
                              "containing configuration data");
DEFINE_string(inter_island_disorder, "none",
              "type of disorder in inter-island interactions: uniform-phase, bernoulli-amplitude"
              ", normal-amplitude");
DEFINE_string(intra_island_disorder, "none",
              "type of disorder in intra-island interactionsi: uniform-phase, bernoulli-amplitude"
              ", normal-amplitude");
DEFINE_bool(proximity_form, false, "coupling interaction proximity form");
DEFINE_double(island_spacing, 0, "island spacing term in proximity form of coupling interactions");
DEFINE_int32(threads, 1, "number of threads for mc sweeps");
DEFINE_string(output_type, "plain", "output format: plain or binary");
DEFINE_int32(print_freq, 1, "print output every print_freq sweeps");
DEFINE_string(array_type, "square", "Lattice array geometry: \"square\" or \"triangular\"");
DEFINE_string(island_type, "square", "Lattice island geometry: \"square\" or \"triangular\"");

// Annealing parameters
DEFINE_string(temp_range, "", "temperature range in the form low,high");
DEFINE_int32(mc_sweeps, 1000, "number of monte-carlo sweeps per annealing step");
DEFINE_int32(annealing_steps, 64, "number of annealing steps from high temperature to low");
DEFINE_int32(overrelax_freq, 10, "How many more over-relaxation steps to perform than metro steps");
DEFINE_double(island_local_field, 0, "Local external field amplitude which each island couples to");
DEFINE_string(island_field_config, "uniform", "Local field configuration for each island: "
                                       "ordered/uniform(parallel to x) or random");
DEFINE_double(spin_local_field, 0, "Local external field amplitude which each spin on each island"
                                   " couples to.");
DEFINE_string(spin_field_config, "uniform", "Local field configuration for each spin: "
                                             "ordered/uniform(parallel to x) or random");

// Hysteresis parameters
DEFINE_int32(steps_per_cycle, 0, "Number of field steps per hysteresis cycle");
DEFINE_int32(hyst_cycles, 1, "Number of hysteresis cycles");
DEFINE_double(temperature, 0, "Fixed temperatures for the hysteresis simulations");

// Output
DEFINE_string(output_file, "", "If set, write output to output_file. "
                               "If not specified, output to stdout.");
DEFINE_bool(overwrite, false, "If set, overwrite the specified output file if it already exists.");

// Visualization parameters
DEFINE_string(window_dims, "512,512", "comma-seperated visual window width,height");
DEFINE_bool(visuals, false, "show visual window of simulations array");

XYChainMPI::Config ParseXYChainConfig(const string& str) {
   if (str == "ordered" || str == "uniform")
      return XYChainMPI::CONFIG_UNIFORM;
   else if (str == "twisted")
      return XYChainMPI::CONFIG_TWISTED;
   else if (str == "twistedpair" || str == "twisted-pair")
      return XYChainMPI::CONFIG_TWISTED_PAIR;
   else 
      return XYChainMPI::CONFIG_RANDOM;
}

XYChainMPI::LatticeType ParseXYChainLatticeType(const string& str) {
   if (str == "square") {
      return XYChainMPI::SQUARE_TYPE;
   } else if (str == "triangular") {
      return XYChainMPI::TRIANGULAR_TYPE;
   } else {
      return XYChainMPI::TRIANGULAR_TYPE;
   }
}

XYChainMPI::Disorder ParseXYChainDisorder(const string& str) {
   if (HasPrefix(str, "uniform-phase"))
      return XYChainMPI::DISORDER_UNIFORM_PHASE;
   else if (HasPrefix(str, "bernoulli-amplitude"))
      return XYChainMPI::DISORDER_BERNOULLI_AMPLITUDE;
   else if (HasPrefix(str, "normal-amplitude"))
      return XYChainMPI::DISORDER_NORMAL_AMPLITUDE;
   else if (HasPrefix(str, "diluted"))
      return XYChainMPI::DISORDER_DILUTED;
   else
      return XYChainMPI::DISORDER_NONE;
}

ChainArgs::OutputType ParseOutputType(const string& str) {
   if (str == "plain" || str == "ascii")
      return ChainArgs::OUTPUT_PLAIN;
   else if (str == "binary" || str == "raw")
      return ChainArgs::OUTPUT_BINARY;
   else
      return ChainArgs::OUTPUT_PLAIN;
}

ChainArgs MakeChainArgs() {
   ChainArgs args;
   
   args.array_dims = ParseIntList(FLAGS_array_dims);
   args.array_type = ParseXYChainLatticeType(FLAGS_array_type);
   args.island_dims = ParseIntList(FLAGS_island_dims);
   args.island_type = ParseXYChainLatticeType(FLAGS_island_type);

   args.n_threads = FLAGS_threads;

   args.initial_config = ParseXYChainConfig(FLAGS_initial_config);

   args.inter_island_disorder = ParseXYChainDisorder(FLAGS_inter_island_disorder);
   args.has_inter_island_disorder_param = ParseDoubleParam(FLAGS_inter_island_disorder,
                                                           &args.inter_island_disorder_param);

   args.intra_island_disorder = ParseXYChainDisorder(FLAGS_intra_island_disorder);
   args.has_intra_island_disorder_param = ParseDoubleParam(FLAGS_intra_island_disorder,
                                                           &args.intra_island_disorder_param);

   args.inter_island_coupling_amp = FLAGS_inter_island_coupling;
   args.intra_island_coupling_amp = FLAGS_intra_island_coupling;

   args.proximity_form = FLAGS_proximity_form;
   args.island_spacing = FLAGS_island_spacing;

   args.print_freq = FLAGS_print_freq;
   args.output_type = ParseOutputType(FLAGS_output_type);

   return args;
}

AnnealArgs MakeAnnealArgs() {
   AnnealArgs args;

   args.chain = MakeChainArgs();

   args.mc_sweeps = FLAGS_mc_sweeps;
   args.annealing_steps = FLAGS_annealing_steps;
   args.overrelax_freq = FLAGS_overrelax_freq;

   vector<double> t_range = ParseDoubleList(FLAGS_temp_range);
   if (t_range.size() != 2) {
      args.temp_low = -1;
      args.temp_high = -1;
   } else {
      args.temp_low = t_range[0];
      args.temp_high = t_range[1];
   }


   args.visuals = FLAGS_visuals;
   vector<int> win_dims = ParseIntList(FLAGS_window_dims);
   if (win_dims.size() != 2) {
      args.window_w = -1;
      args.window_h = -1;
   } else {
      args.window_w = win_dims[0];
      args.window_h = win_dims[1];
   }

   args.arr_field = FLAGS_island_local_field;
   args.arr_field_type = ParseXYChainConfig(FLAGS_island_field_config);
   args.isl_field = FLAGS_spin_local_field;
   args.isl_field_type = ParseXYChainConfig(FLAGS_spin_field_config);

   return args;
}

HystArgs MakeHystArgs() {
   HystArgs args;

   args.chain = MakeChainArgs();

   args.mc_sweeps = FLAGS_mc_sweeps;

   args.field_amp = FLAGS_spin_local_field;
   args.field_steps = FLAGS_steps_per_cycle;
   args.field_cycles = FLAGS_hyst_cycles;

   args.temp = FLAGS_temperature;

   args.visuals = FLAGS_visuals;
   vector<int> win_dims = ParseIntList(FLAGS_window_dims);
   if (win_dims.size() != 2) {
      args.window_w = -1;
      args.window_h = -1;
   } else {
      args.window_w = win_dims[0];
      args.window_h = win_dims[1];
   }

   return args;
}

PTArgs MakePTArgs() {
   PTArgs args;

   // TODO(jw): implement me

   return args;
}

class AnnealingCommand : public Command {
 public:
   AnnealingCommand() {
   }

   int Run() {
      AnnealArgs args = MakeAnnealArgs();
      if (!VerifyAnnealArgs(args)) {
         fprintf(stderr, "Invalid args\n");
         return 1;
      }
      return RunAnnealing(args);
   }

   int RunHelp() {
      fprintf(stderr, "xychain annealing [flags...]\n"
            "Collect observable values versus temperature. The temperature is initially relatively\n"
            "high, and gradually lowered to some relatively low value, performing monte-carlo sweeps\n"
            "at each step along the way.\n\n"
            "Flags:\n"
            "  --temp_range: a comma-seperated pair of temperatures. Simulation will\n"
            "        start at temp_high, and gradually the temperature is lowered to temp_low.\n"
            "  --mc_sweeps: number of monte-carlo sweeps to perform per temperature step.\n"
            "  --annealing_steps: number of temperature steps between temp_low and temp_high.\n"
            "  --overrelax_freq: only perform a metropolis sweep once every overrelax_freq steps.\n"
            "        Higher numbers lead to a greater ratio of over-relaxation sweeps.\n"
            "        Setting overrelax_freq=1 disables over relaxation sweeps.\n"
            "  --island_local_field: amplitude for the local external field each island couples to.\n"
            "  --island_field_config: configuration for local field of each island\n"
            "        (uniform or random)\n"
            "  --spin_local_field: amplitude for the local external field each spin couples to.\n"
            "  --spin_field_config: configuration for the local field of each spin\n"
            "        (uniform or random)\n"
            "Additionally, all the flags which describe the array of lattices apply. These\n"
            "are listed by the \"xychain help\" general usage.\n\n");
      return 1;
   }
};

class HysteresisCommand : public Command {
 public:
   HysteresisCommand() {
   }

   int Run() {
      HystArgs args = MakeHystArgs();
      if (!VerifyHystArgs(args)) {
         fprintf(stderr, "Invalid args\n");
         return 1;
      }
      return RunHysteresis(args);
   }

   int RunHelp() {
      fprintf(stderr, "xychain hyst [flags...]\n"
            "Run a hysteresis simulation, in which temperature is held\n"
            "constant and the external field is varied as H * cos(t)\n\n"
            "Flags:\n"
            "  --steps_per_cycle: the number of times the external field is\n"
            "     varied per hysteresis cycle\n"
            "  --hyst_cycles: the number of hysteresis cycles to perform\n"
            "  --temperature: temperature to simulate\n"
            "  --mc_sweeps: number of metropolis sweeps to perform per step\n"
            "  --spin_local_field: field amplitude H\n\n");
      return 1;
   }
};

class PTemperingCommand : public Command {
 public:
   PTemperingCommand() {
   }

   int Run() {
      PTArgs args = MakePTArgs();
      if (!VerifyPTArgs(args)) {
         fprintf(stderr, "Invalid args\n");
         return 1;
      }
      return RunPTempering(args);
   }

   int RunHelp() {
      fprintf(stderr, "xychain ptempering [flags...]\n"
            "TODO(jw): Impliment me\n");
      return 1;
   }
};

map<string, Command*> xychain_commands;
void InitCommands() {
   xychain_commands["annealing"] = (Command*)(new AnnealingCommand());
//   xychain_commands["ptempering"] = (Command*)(new PTemperingCommand());
   xychain_commands["hyst"] = (Command*)(new HysteresisCommand());
}

int PrintUsage() {
   fprintf(stderr, "XY Chain model simulator.\n"
         "Simulates a rectangular array of rectangular islands, \n"
         "where each island is a typical XY lattice with \n"
         "O2 vector spins coupled to nearest-neighbors.\n"
         "Additionally, the mean-field of islands are also \n"
         "coupled to nearest-neighbor islands.\n\n"
         "Must be run under MPI to function properly. Additionally,\n"
         "the number of processes must be at least 2.\n"
         "\n"
         "Valid commands are:\n");
   for (const pair<string, Command*>& kv : xychain_commands) {
      fprintf(stderr, "  %s\n", kv.first.c_str());
   }
   fprintf(stderr, "  help\n"
         "\n"
         "To see basic usage for a command, run \"xychain help [command]\"\n\n");
   fprintf(stderr, "Common flags for specifying chain/array and lattice parameters:\n"
         "  --flagfile: a file to read flag values from.\n"
         "  --array_dims: comma-seperated list of array dimensions.\n"
         "        Should have length >= 1.\n"
         "  --array_type: lattice array geometry, either \"square\" or \"triangular\"\n"
         "  --island_dims: comma-seperated list of dimensions for each island in the array.\n"
         "        Should have length >= 1.\n"
         "  --island_type: lattice island geometry, either \"square\" or \"triangular\"\n"
         "  --inter_island_coupling: the coupling amplitude between islands.\n"
         "  --intra_island_coupling: the coupling amplitude between spins on an island.\n"
         "  --initial_config: type of initial configuration for spins\n"
         "        ex: \"uniform\" or \"random\".\n"
         "  --inter_island_disorder: type of (random) disorder in inter-island interactions\n"
         "        \"uniform-phase\", \"bernoulli-amplitude\", \"normal-amplitude\"\n"
         "        \"diluted[p]\"\n"
         "  --intra_island_disorder: type of (random) disorder in intra-island interactions\n"
         "        as with inter_island_disorder\n"
         "  --proximity_form: whether or not to use proximity form coupling amplitudes\n"
         "  --island_spacing: if using proximity_form couplings, the island spacing parameter\n"
         "        for inter-island interactions.\n"
         "  --threads: number of threads to use for sweeps. Note that the array dimensions\n"
         "        are divided among threads, so small arrays may not be fully parallelizable.\n"
         "  --print_freq: print output only every print_freq sweeps.\n"
         "  --output_type: \"ascii\" or \"binary\" output format. Binary format has higher\n"
         "        decimal precision.\n\n");
   return 1;
}

int main(int argc, char** argv) {
   MPI::Init(argc, argv);

   google::SetUsageMessage("Usage: " + string(argv[0]) + " [...]");
   google::ParseCommandLineFlags(&argc, &argv, true);

   InitCommands();

   if (argc < 2) {
      return PrintUsage();
   }

   string command(argv[1]);
   // Run the requested command
   if (xychain_commands.find(command) != xychain_commands.end()) {
      return xychain_commands[command]->Run();
   // Run the "help" command
   } else if (command == "help") {
      if (argc < 3) {
         return PrintUsage();
      } else {
         // Run the "help (command)" command
         string query(argv[2]);
         if (xychain_commands.find(query) != xychain_commands.end()) {
            return xychain_commands[query]->RunHelp();
         } else {
            return PrintUsage();
         }
      }
   }

   fprintf(stderr, "Unknown command: %s", command.c_str());
   return PrintUsage();
}
