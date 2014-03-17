#ifndef XYCHAIN_H
#define XYCHAIN_H

#include "rng.h"
#include "spins.h"
#include "threadpool.h"
#include <vector>

class XYChainMPI {
 public:
   enum MPIMessageTags {
      SLICE_TAG,
      SLICE_BUFLEN_TAG,
      SLICE_META_TAG,
      SLICE_META_BUFLEN_TAG,
      REPORT_TAG,
      REPORTS_DONE_TAG,
      COMMAND_TAG
   };
   enum RemoteCommands {
      COMMAND_DERP,
      COMMAND_EXIT,
      COMMAND_METRO_SWEEP,
      COMMAND_OVERRELAX_SWEEP,
      COMMAND_RESYNC,
      COMMAND_TEST
   };
   struct Measures {
      double energy;
      spins::O2 total_spin;
 //     std::vector<double> spin_cos;
//      std::vector<double> spin_sin;
      double m_isl_avg;
      double m2_isl_avg;
   };
   enum Config {
      CONFIG_RANDOM,
      CONFIG_UNIFORM,
      CONFIG_TWISTED,
      CONFIG_TWISTED_PAIR
//      CONFIG_SERIAL_IN
   };
   enum Disorder {
      DISORDER_NONE,
      DISORDER_BERNOULLI_AMPLITUDE,
      DISORDER_NORMAL_AMPLITUDE,
      DISORDER_UNIFORM_PHASE,
      DISORDER_DILUTED
   };
   enum BoundaryCond {
      PERIODIC_BOUNDS
   };
   enum LatticeType {
      SQUARE_TYPE,
      TRIANGULAR_TYPE
   };

   XYChainMPI(const std::vector<int>& array_dimensions,
              const std::vector<int>& island_dimensions,
              const LatticeType& array_type,
              const LatticeType& island_type);

   ~XYChainMPI();

   static int GetNumNeighbors(const std::vector<int>& dims, const LatticeType& type);

   void SetBeta(double beta) {fBeta = beta;}
   double GetBeta() const {return fBeta;}
   int GetNarr() const {return fNarr;}
   int GetIslandVolume() const {return fIslVol;}
   int GetTotalVolume() const {return fTotalVol;}

   void MetroSweep();
   void OverRelaxSweep();

   void Init(XYChainMPI::Config initial_config);
   void InitArrayCoupling(double r, Disorder type);
   void InitArrayCoupling(double r, Disorder type, double disorder_param);
   void InitIslandCoupling(double r, Disorder type);
   void InitIslandCoupling(double r, Disorder type, double disorder_param);
   void ScaleArrayCouplingAmplitude(double r);
   void ScaleIslandCouplingAmplitude(double r);

   void InitIslandField(double r, Config type);
   void InitArrayField(double r, Config type);

   XYChainMPI::Measures GetMeasures() {
      if (fIsDirty) {
         fMeasures = GetMeasuresExhaustive();
         fIsDirty = false;
      }
      return fMeasures;
   }
   XYChainMPI::Measures GetIslandMeasures(int i) {
      if (fIsDirty) {
         fMeasures = GetMeasuresExhaustive();
         fIsDirty = false;
      }
      return fIslands[i].measures;
   }
   XYChainMPI::Measures GetMeasuresExhaustive() const;

   static void AddMeasures(const XYChainMPI::Measures& m1,
                           XYChainMPI::Measures* m2);
   static void ZeroMeasures(XYChainMPI::Measures* m);

   void FlattenStaticInfo(int arr_i_begin, int arr_i_end, char* buf);
   void UnflattenStaticInfo(const char* buf);
   void ReceiveStaticInfo(int src);
   void SendStaticInfo(int arr_i_begin, int arr_i_end, int dest);

   void FlattenSlice(int arr_i_begin, int arr_i_end, char* buf);
   void UnflattenSlice(const char* buf);
   void ReceiveSlice(int src);
   void SendSlice(int arr_i_begin, int arr_i_end, int dest);

   void SendMeasureDeltas(const XYChainMPI::Measures* meas);

   void MetroSweepMPI();
   void OverRelaxSweepMPI();
   void ResyncMPI();
   void CommandLoop();

   static XYChainMPI::Measures MasterMetroSweeps(int n_sweeps);
   static XYChainMPI::Measures MasterOverRelaxSweeps(int n_sweeps);
   static XYChainMPI::Measures MasterResync();
   static XYChainMPI::Measures MasterCollect(int n_sweeps);
   static void MasterShutdownSlaves();
   static void MasterPrintTest();

   // Testing
   int GetArrayNeighborIndex(int index, int neighbor_index) {
      return fIslands[index].neighbor_ind[neighbor_index];
   }
   spins::O2Coupling GetArrayNeighborCoupling(int index, int neighbor_index) {
      return fIslands[index].neighbor_coupling[neighbor_index];
   }
   int GetIslandSiteNeighborIndex(int index, int isl_index, int neighbor_index) {
      return fIslands[index].sites[isl_index].neighbor_ind[neighbor_index];
   }
   spins::O2Coupling GetIslandSiteNeighborCoupling(int index, int isl_index, int neighbor_index) {
      return fIslands[index].sites[isl_index].neighbor_coupling[neighbor_index];
   }

//   static XYChain* Clone(const XYChain& other);

   #ifdef VISUALS
   void GLRender(double x, double y, double scale);
   void GLRender(double scale);
   #endif // VISUALS

 private:
   void InitThreading();
   void InitIslands(const LatticeType& arr_type, const LatticeType& isl_type);
   void ClearThreading();
   void ClearIslands();

   struct SpinSite {
      spins::O2 spin;
      spins::O2Coupling* neighbor_coupling;
      int* neighbor_ind;
      spins::O2 local_field;
   };

   struct Island {
      XYChainMPI::SpinSite* sites;
      spins::O2Coupling* neighbor_coupling;
      int* neighbor_ind;
      XYChainMPI::Measures measures;
      spins::O2 local_field;
   };

   struct ThreadArgs {
      XYChainMPI* self;
      int slice_begin;
      int slice_end;
      rng::RandomGen rng;
      int n_accept;
      int n_reject;
      double window;
      XYChainMPI::Measures measures_ret;
   };

   struct sum_meas_change_args {
      XYChainMPI* self;
      ThreadPool* pool;
      int slice_0;
      int slice_1;
      XYChainMPI::Measures* change;
   };

   void SweepMPI(void (*sweep_func)(void* arg));
   bool DoRemoteCommand();
   static void MasterSendCommand(int cmd, int num);
   XYChainMPI::Measures GetMeasuresExhSlice(int arr_i_begin, int arr_i_end) const;

   void Sweep(void (*sweep_func)(void* arg));
   void SumMeasureChanges();
   static void MetroSweepSlice(void* args);
   static void OverRelaxSweepSlice(void* args);
   static void SumMeasureChangesRec(void* args);

   static int* MakeNeighborIndices(const std::vector<int>& dims, int index,
                                   const LatticeType& type);
   static int* MakeNeighborIndicesSquare(const std::vector<int>& dims, int index);
   static int* MakeNeighborIndicesTriangle(const std::vector<int>& dims, int index);

public:
   static spins::O2 RandomUnitO2(rng::RandomGen* rng);
   static spins::O2Coupling RandomCoupling(rng::RandomGen* rng, double r, Disorder type);
   static spins::O2Coupling RandomCoupling(rng::RandomGen* rng, double r, Disorder type,
                                           double param);
private:
   const int fNthreads;
//   ThreadPool fJobPool;
   XYChainMPI::ThreadArgs* fThreadArgs;
   int fNslice;

   std::vector<int> fArrayDims;
   std::vector<int> fIslandDims;
   // Should really be const
   int fNarr;
   int fIslVol;
   int fTotalVol;
   XYChainMPI::Island* fIslands;

   XYChainMPI::Measures fMeasures;
   double fBeta;

   int fNaccept;
   int fNreject;
   double fWindow;

   // If true, we need to recalculate measures
   // exhaustively.
   bool fIsDirty;

   const int fNMpiIds;
   char* fSendBuf;

   const int fArrayCoord;
   const int fIslandCoord;
};

#endif
