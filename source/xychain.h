#ifndef XYCHAIN_H
#define XYCHAIN_H

#include "rng.h"
#include "spins.h"
#include "threadpool.h"
#include <vector>

class XYChain {
 public:
   struct Measures {
      double energy;
      spins::O2 total_spin;
      std::vector<double> spin_cos;
      std::vector<double> spin_sin;
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

   XYChain(std::vector<int> array_dimensions,
           std::vector<int> island_dimensions,
           int n_threads);

   ~XYChain();

   void SetBeta(double beta) {fBeta = beta;}
   double GetBeta() const {return fBeta;}
   int GetNarr() const {return fNarr;}
   int GetIslandVolume() const {return fIslVol;}
   int GetTotalVolume() const {return fTotalVol;}

   void MetroSweep();
   void OverRelaxSweep();

   void Init(XYChain::Config initial_config);
   void InitArrayCoupling(double r, Disorder type);
   void InitArrayCoupling(double r, Disorder type, double disorder_param);
   void InitIslandCoupling(double r, Disorder type);
   void InitIslandCoupling(double r, Disorder type, double disorder_param);
   void ScaleArrayCouplingAmplitude(double r);
   void ScaleIslandCouplingAmplitude(double r);

   void InitIslandField(double r, Config type);
   void InitArrayField(double r, Config type);

   XYChain::Measures GetMeasures() {
      if (fIsDirty) {
         fMeasures = GetMeasuresExhaustive();
         fIsDirty = false;
      }
      return fMeasures;
   }
   XYChain::Measures GetIslandMeasures(int i) {
      if (fIsDirty) {
         fMeasures = GetMeasuresExhaustive();
         fIsDirty = false;
      }
      return fIslands[i].measures;
   }
   XYChain::Measures GetMeasuresExhaustive() const;

   static void AddMeasures(const XYChain::Measures& m1,
                           XYChain::Measures* m2);
   static void ZeroMeasures(XYChain::Measures* m);

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
   void InitIslands();
   void ClearThreading();
   void ClearIslands();

   struct SpinSite {
      spins::O2 spin;
      spins::O2Coupling* neighbor_coupling;
      int* neighbor_ind;
      spins::O2 local_field;
   };

   struct Island {
      XYChain::SpinSite* sites;
      spins::O2Coupling* neighbor_coupling;
      int* neighbor_ind;
      XYChain::Measures measures;
      spins::O2 local_field;
   };

   struct ThreadArgs {
      XYChain* self;
      int slice_begin;
      int slice_end;
      rng::RandomGen rng;
      int n_accept;
      int n_reject;
      double window;
      XYChain::Measures measures_ret;
   };

   struct sum_meas_change_args {
      XYChain* self;
      ThreadPool* pool;
      int slice_0;
      int slice_1;
      XYChain::Measures* change;
   };

   void Sweep(void (*sweep_func)(void* arg));
   void SumMeasureChanges();
   static void MetroSweepSlice(void* args);
   static void OverRelaxSweepSlice(void* args);
   static void SumMeasureChangesRec(void* args);

   static int* MakeNeighborIndices(const std::vector<int>& dims, int index);

public:
   static spins::O2 RandomUnitO2(rng::RandomGen* rng);
   static spins::O2Coupling RandomCoupling(rng::RandomGen* rng, double r, Disorder type);
   static spins::O2Coupling RandomCoupling(rng::RandomGen* rng, double r, Disorder type,
                                           double param);
private:
   const int fNthreads;
   ThreadPool fJobPool;
   XYChain::ThreadArgs* fThreadArgs;
   int fNslice;

   std::vector<int> fArrayDims;
   std::vector<int> fIslandDims;
   // Should really be const
   int fNarr;
   int fIslVol;
   int fTotalVol;
   XYChain::Island* fIslands;

   XYChain::Measures fMeasures;
   double fBeta;

   int fNaccept;
   int fNreject;
   double fWindow;

   // If true, we need to recalculate measures
   // exhaustively.
   bool fIsDirty;
};

#endif
