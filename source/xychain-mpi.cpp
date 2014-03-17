#include "xychain-mpi.h"
#include <math.h>
#include <unordered_set>
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;
using std::min;

#include <mpi.h>

#include "notification.h"

#ifdef VISUALS
 #ifdef LINUX
  #include <GL/gl.h>
 #endif  // LINUX

 #ifdef OSX
  #include <OpenGL/gl.h>
 #endif  // OSX
#endif  // VISUALS

int XYChainMPI::GetNumNeighbors(const vector<int>& dims, const LatticeType& type) {
   if (type == SQUARE_TYPE) {
      return 2*dims.size();
   } else if (type == TRIANGULAR_TYPE) {
      int pow = 1;
      for (int i = 0; i < dims.size(); ++i) {
         pow *= 2;
      }
      return 2*(pow - 1);
   }
   return 0;
}

XYChainMPI::XYChainMPI(const vector<int>& array_dimensions,
                       const vector<int>& island_dimensions,
                       const LatticeType& array_type,
                       const LatticeType& island_type)
      : fNthreads(MPI::COMM_WORLD.Get_size()-1),
//        fJobPool(fNthreads),
        fArrayDims(array_dimensions),
        fIslandDims(island_dimensions),
        fBeta(1),
        fNaccept(0),
        fNreject(0),
        fWindow(1),
        fNMpiIds(MPI::COMM_WORLD.Get_size()),
        fArrayCoord(GetNumNeighbors(fArrayDims, array_type)),
        fIslandCoord(GetNumNeighbors(fIslandDims, island_type)) {
   InitIslands(array_type, island_type);
   InitThreading();
   fIsDirty = true;

   const int arr_coord_num = fArrayCoord;
   const int isl_coord_num = fIslandCoord;
   const int n_arr_slice = fNarr / fArrayDims.back();
   const int n_spins = n_arr_slice * fIslVol;
   const int maxbuflen = 1 + 2 * sizeof(int)
                        + n_arr_slice * sizeof(spins::O2)
                        + n_arr_slice * arr_coord_num * sizeof(spins::O2Coupling)
                        + n_arr_slice * arr_coord_num * sizeof(int)
                        + n_spins * sizeof(spins::O2)
                        + n_spins * isl_coord_num * sizeof(spins::O2Coupling)
                        + n_spins * isl_coord_num * sizeof(int);
   fSendBuf = new char[maxbuflen];
}

void XYChainMPI::InitIslands(const LatticeType& arr_type, const LatticeType& isl_type) {
   fNarr = 1;
   for (int i = 0; i < fArrayDims.size(); ++i) {
      fNarr *= fArrayDims[i];
   }

   fIslVol = 1;
   for (int i = 0; i < fIslandDims.size(); ++i) {
      fIslVol *= fIslandDims[i];
   }

   fTotalVol = fNarr * fIslVol;

   fIslands = new XYChainMPI::Island[fNarr];
   int arr_coord_num = fArrayCoord;
   int isl_coord_num = fIslandCoord;
   for (int i = 0; i < fNarr; ++i) {
      fIslands[i].sites = new XYChainMPI::SpinSite[fIslVol];
      for (int j = 0; j < fIslVol; ++j) {
         fIslands[i].sites[j].neighbor_coupling = new spins::O2Coupling[isl_coord_num];
         fIslands[i].sites[j].neighbor_ind = MakeNeighborIndices(fIslandDims, j, isl_type);
         fIslands[i].sites[j].local_field = {0, 0};
      }
      fIslands[i].neighbor_coupling = new spins::O2Coupling[arr_coord_num];
      fIslands[i].neighbor_ind = MakeNeighborIndices(fArrayDims, i, arr_type);
      fIslands[i].local_field = {0, 0};
   }

   fIsDirty = true;
}

void XYChainMPI::InitThreading() {
   fNslice = fNthreads * 2;
   if (fArrayDims.back() < fNslice) {
      fNslice = fArrayDims.back();
   }
   fThreadArgs = new XYChainMPI::ThreadArgs[fNslice];
   double slice_size = (double)fArrayDims.back() / fNslice;
   for (int i = 0; i < fNslice; ++i) {
      fThreadArgs[i].self = this;
      fThreadArgs[i].slice_begin = (int)(i * slice_size);
      fThreadArgs[i].slice_end = (int)((i+1) * slice_size);
      fThreadArgs[i].n_accept = 0;
      fThreadArgs[i].n_reject = 0;
      fThreadArgs[i].window = 1;
      ZeroMeasures(&fThreadArgs[i].measures_ret);
   }
}

void XYChainMPI::ClearThreading() {
   delete [] fThreadArgs;
}

void XYChainMPI::ClearIslands() {
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         delete [] fIslands[i].sites[j].neighbor_coupling;
         delete [] fIslands[i].sites[j].neighbor_ind;
      }
      delete [] fIslands[i].sites;
      delete [] fIslands[i].neighbor_coupling;
      delete [] fIslands[i].neighbor_ind;
   }
   delete [] fIslands;
}

XYChainMPI::~XYChainMPI() {
   ClearIslands();
   ClearThreading();

   delete [] fSendBuf;
}

int* XYChainMPI::MakeNeighborIndices(const std::vector<int>& dims, int index, const LatticeType& type) {
   if (type == SQUARE_TYPE) {
      return MakeNeighborIndicesSquare(dims, index);
   } else if (type == TRIANGULAR_TYPE) {
      return MakeNeighborIndicesTriangle(dims, index);
   } else {
      return 0;
   }
}

int CoordsToIndex(const vector<int>& dims, const vector<int>& coords) {
   int idx = 0;
   int pow = 1;
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      idx += pow * coords[dim_i];
      pow *= dims[dim_i];
   }
   return idx;
}

int* XYChainMPI::MakeNeighborIndicesTriangle(const std::vector<int>& dims, int index) {
   int coord_num = GetNumNeighbors(dims, TRIANGULAR_TYPE);
   int* neighbor_i = new int[coord_num];
   vector<int> index_coords(dims.size());
   int ind = index;
   // Get the coordinates of the given index
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      index_coords[dim_i] = ind % dims[dim_i];
      ind /= dims[dim_i];
   }
   
   // Need to iterate over all possible combinations of dimensions...
   for (int bitmask = 0; bitmask < coord_num/2; ++bitmask) {
      vector<int> neighbor_coords = index_coords;
      for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
         // If bit dim_i is set in bitmask, then add 1 to coords[dim_i] mod periodic boundary conditions
         if ((bitmask+1) & (1 << dim_i)) {
            neighbor_coords[dim_i] = (neighbor_coords[dim_i]+1) % dims[dim_i];
         }
      }
      neighbor_i[bitmask] = CoordsToIndex(dims, neighbor_coords);
   }

   // Now in reverse direction
   for (int bitmask = 0; bitmask < coord_num/2; ++bitmask) {
      vector<int> neighbor_coords = index_coords;
      for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
         if ((bitmask+1) & (1 << dim_i)) {
            neighbor_coords[dim_i] = (neighbor_coords[dim_i] + dims[dim_i] - 1) % dims[dim_i];
         }
      }
      neighbor_i[bitmask + coord_num/2] = CoordsToIndex(dims, neighbor_coords);
   }

   return neighbor_i;
}

int* XYChainMPI::MakeNeighborIndicesSquare(const std::vector<int>& dims, int index) {
   int coord_num = GetNumNeighbors(dims, SQUARE_TYPE);
   int* neighbor_i = new int[coord_num];
   vector<int> index_coords(dims.size());
   int ind = index;
   // Get the coordinates of the given index
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      index_coords[dim_i] = ind % dims[dim_i];
      ind /= dims[dim_i];
   }
   // Get the indices of the neighboring coordinates
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      index_coords[dim_i] = (index_coords[dim_i]+1) % dims[dim_i];
      neighbor_i[dim_i] = CoordsToIndex(dims, index_coords);
      index_coords[dim_i] = (index_coords[dim_i]+dims[dim_i]-1) % dims[dim_i];
   }
   // Do it again, for negative neighbors
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      index_coords[dim_i] = (index_coords[dim_i]+dims[dim_i]-1) % dims[dim_i];
      neighbor_i[dim_i+dims.size()] = CoordsToIndex(dims, index_coords);
      index_coords[dim_i] = (index_coords[dim_i]+1) % dims[dim_i];
   }
   // We're done!
   return neighbor_i;
}

spins::O2 XYChainMPI::RandomUnitO2(rng::RandomGen* rng) {
   double pi = acos(0)*2;
   double angle = rng->RandomUniform()*2*pi;
   double x = cos(angle);
   double y = sin(angle);
   return {x,y};
}

void XYChainMPI::Init(XYChainMPI::Config initial_config) {
   rng::RandomGen rng;
   int arr_coord_num = fArrayCoord;
   for (int i = 0; i < fNarr; ++i) {
      // Initialize constant coupling (identity matrix)
      for (int j = 0; j < arr_coord_num; ++j) {
         fIslands[i].neighbor_coupling[j] = spins::O2Coupling_I;
      }

      spins::O2 total_spin = spins::O2_zero;
      // Now initialize each lattice site in the island
      for (int site_i = 0; site_i < fIslVol; ++site_i) {
         spins::O2 spinval;
         if (initial_config == CONFIG_RANDOM) {
            spinval = RandomUnitO2(&rng);
         } else if (initial_config == CONFIG_UNIFORM) {
            spinval = spins::O2_i;
//         } else if (initial_config == CONFIG_SERIAL_IN) {
//            spinval = fIslands[i].sites[site_i].spin;
         } else if (initial_config == CONFIG_TWISTED) {
            int x = i % fArrayDims[0];
            double pi = acos(0)*2;
            double angle = 2*pi / fArrayDims[0] * x;
            spinval = {cos(angle), sin(angle)};
         } else if (initial_config == CONFIG_TWISTED_PAIR) {
            int x = i % fArrayDims[0];
            double pi = acos(0)*2;
            double angle = 2*pi / fArrayDims[0] * 2 * x;
            if (x > fArrayDims[0]/2) angle *= -1;
            spinval = {cos(angle), sin(angle)};
         } else {
            spinval = RandomUnitO2(&rng);
         }
         fIslands[i].sites[site_i].spin = spinval;
         spins::O2SumAccum(spinval, &total_spin);

         // Initialize spin coupling on the island to uniform (identity matrix)
         int isl_coord_num = fIslandCoord;
         for (int coupl_i = 0; coupl_i < isl_coord_num; ++coupl_i) {
            fIslands[i].sites[site_i].neighbor_coupling[coupl_i] =
               spins::O2Coupling_I;
         }
      }
      fIslands[i].measures.total_spin = total_spin;
   }

   fIsDirty = true;
}

void XYChainMPI::ScaleArrayCouplingAmplitude(double r) {
   int arr_coord_num = fArrayCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < arr_coord_num; ++j) {
         spins::O2ScaleCouplingAccum(r, &fIslands[i].neighbor_coupling[j]);
      }
   }

   fIsDirty = true;
}

void XYChainMPI::ScaleIslandCouplingAmplitude(double r) {
   int isl_coord_num = fIslandCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         for (int k = 0; k < isl_coord_num; ++k) {
            spins::O2ScaleCouplingAccum(r, 
                  &fIslands[i].sites[j].neighbor_coupling[k]);
         }
      }
   }

   fIsDirty = true;
}

void XYChainMPI::InitArrayCoupling(double r, XYChainMPI::Disorder disorder) {
   rng::RandomGen rng;
   int arr_coord_num = fArrayCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int  j = 0; j < arr_coord_num/2; ++j) {
         spins::O2Coupling c = RandomCoupling(&rng, r, disorder);
         fIslands[i].neighbor_coupling[j] = c;
         int ni = fIslands[i].neighbor_ind[j];
         fIslands[ni].neighbor_coupling[j+arr_coord_num/2] = spins::O2CouplingTranspose(c);
      }
   }

   fIsDirty = true;
}

void XYChainMPI::InitArrayCoupling(double r, XYChainMPI::Disorder disorder, double param) {
   rng::RandomGen rng;
   int arr_coord_num = fArrayCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int  j = 0; j < arr_coord_num/2; ++j) {
         spins::O2Coupling c = RandomCoupling(&rng, r, disorder, param);
         fIslands[i].neighbor_coupling[j] = c;
         int ni = fIslands[i].neighbor_ind[j];
         fIslands[ni].neighbor_coupling[j+arr_coord_num/2] = spins::O2CouplingTranspose(c);
      }
   }

   fIsDirty = true;
}

void XYChainMPI::InitIslandField(double r, XYChainMPI::Config type) {
   rng::RandomGen rng;
   int isl_coord_num = fIslandCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         spins::O2 f = spins::O2_i;
         if (type == XYChainMPI::CONFIG_RANDOM) {
            f = RandomUnitO2(&rng);
         }
         spins::O2ScaleAccum(r, &f);
         fIslands[i].sites[j].local_field = f;
      }
   }

   fIsDirty = true;
}

void XYChainMPI::InitArrayField(double r, XYChainMPI::Config type) {
   rng::RandomGen rng;
   int arr_coord_num = fArrayCoord;
   for (int i = 0; i < fNarr; ++i) {
      spins::O2 f = spins::O2_i;
      if (type == XYChainMPI::CONFIG_RANDOM) {
         f = RandomUnitO2(&rng);
      }
      spins::O2ScaleAccum(r, &f);
      fIslands[i].local_field = f;
   }

   fIsDirty = true;
}

void XYChainMPI::InitIslandCoupling(double r, XYChainMPI::Disorder disorder) {
   rng::RandomGen rng;
   int isl_coord_num = fIslandCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         for (int k = 0; k < isl_coord_num/2; ++k) {
            spins::O2Coupling c = RandomCoupling(&rng, r, disorder);
            fIslands[i].sites[j].neighbor_coupling[k] = c;
            int nj = fIslands[i].sites[j].neighbor_ind[k];
            fIslands[i].sites[nj].neighbor_coupling[k+isl_coord_num/2] = spins::O2CouplingTranspose(c);
         }
      }
   }

   fIsDirty = true;
}

void XYChainMPI::InitIslandCoupling(double r, XYChainMPI::Disorder disorder, double param) {
   rng::RandomGen rng;
   int isl_coord_num = fIslandCoord;
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         for (int k = 0; k < isl_coord_num/2; ++k) {
            spins::O2Coupling c = RandomCoupling(&rng, r, disorder, param);
            fIslands[i].sites[j].neighbor_coupling[k] = c;
            int nj = fIslands[i].sites[j].neighbor_ind[k];
            fIslands[i].sites[nj].neighbor_coupling[k+isl_coord_num/2] = spins::O2CouplingTranspose(c);
         }
      }
   }

   fIsDirty = true;
}

spins::O2Coupling XYChainMPI::RandomCoupling(rng::RandomGen* rng, double r, Disorder type) {
   switch (type) {
      case DISORDER_NONE:
         return spins::O2ScaleCoupling(r, spins::O2Coupling_I);
      case DISORDER_BERNOULLI_AMPLITUDE:
         return spins::O2ScaleCoupling((rng->RandomUniform() < 0.5 ? r : -r),
                                       spins::O2Coupling_I);
      case DISORDER_NORMAL_AMPLITUDE:
         return spins::O2ScaleCoupling(r * rng->RandomNormal(), spins::O2Coupling_I);
      case DISORDER_UNIFORM_PHASE:
         {
            double pi = acos(0)*2;
            double angle = rng->RandomUniform()*2*pi;
            spins::O2Coupling c = {cos(angle), -sin(angle), sin(angle), cos(angle)};
            return spins::O2ScaleCoupling(r, c);
         }
   }
}

spins::O2Coupling XYChainMPI::RandomCoupling(rng::RandomGen* rng, double r, Disorder type,
                                          double param) {
   switch (type) {
      case DISORDER_NONE:
         return spins::O2ScaleCoupling(r, spins::O2Coupling_I);
      case DISORDER_BERNOULLI_AMPLITUDE:
         return spins::O2ScaleCoupling((rng->RandomUniform() < param ? r : -r),
                                       spins::O2Coupling_I);
      case DISORDER_NORMAL_AMPLITUDE:
         return spins::O2ScaleCoupling(param + r * rng->RandomNormal(), spins::O2Coupling_I);
      case DISORDER_UNIFORM_PHASE:
         {
            double pi = acos(0)*2;
            double angle = (2*rng->RandomUniform()-1)*pi*param;
            spins::O2Coupling c = {cos(angle), -sin(angle), sin(angle), cos(angle)};
            return spins::O2ScaleCoupling(r, c);
         }
      case DISORDER_DILUTED:
         return spins::O2ScaleCoupling((rng->RandomUniform() < param ? r : 0),
                                    spins::O2Coupling_I);
   }
}

XYChainMPI::Measures XYChainMPI::GetMeasuresExhSlice(int arr_i_begin, int arr_i_end) const {
   double total_energy = 0;
   spins::O2 total_spin = spins::O2_zero;
   vector<double> total_cos(fArrayDims.size());
   vector<double> total_sin(fArrayDims.size());
   double m_isl_avg = 0;
   double m2_isl_avg = 0;

   int arr_coord_num = fArrayCoord;
   int isl_coord_num = fIslandCoord;

   // Accumulate island total spin
   for (int i = arr_i_begin; i <= arr_i_end; ++i) {
      int ind = i % fNarr;
      spins::O2 isl_total_spin = spins::O2_zero;
      for (int j = 0; j < fIslVol; ++j) {
         const spins::O2 s = fIslands[ind].sites[j].spin;
         spins::O2SumAccum(s, &isl_total_spin);
      }
      fIslands[ind].measures.total_spin = isl_total_spin;
      double m2 = spins::O2Product(isl_total_spin, isl_total_spin);
      spins::O2SumAccum(isl_total_spin, &total_spin);
      m_isl_avg += sqrt(m2) / fNarr;
      m2_isl_avg += m2 / fNarr;
   }

   // Accumulate inter-island coupling energy term
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      const spins::O2 isl_total_spin = fIslands[i].measures.total_spin;
      for (int arr_adj = 0; arr_adj < arr_coord_num/2; ++arr_adj) {
         const int neighbor_i = fIslands[i].neighbor_ind[arr_adj];
         const spins::O2 isl_adj_total_spin = fIslands[neighbor_i].measures.total_spin;
         const spins::O2Coupling isl_coupling = fIslands[i].neighbor_coupling[arr_adj];
         total_energy -= spins::O2Product(isl_total_spin, isl_coupling, isl_adj_total_spin);
         total_energy -= spins::O2Product(isl_total_spin, fIslands[i].local_field);

         const spins::O2 isl_adj_rot = spins::O2Product(
               spins::O2Coupling_rotate90, isl_adj_total_spin);
//         total_cos[arr_adj] += spins::O2Product(
//               isl_total_spin, isl_coupling, isl_adj_total_spin);
//         total_sin[arr_adj] += spins::O2Product(
//               isl_total_spin, isl_coupling, isl_adj_rot);
      }
   }

   // Accumulate intra-island coupling energy terms
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         const spins::O2 s_a = fIslands[i].sites[j].spin;
         for (int isl_adj = 0; isl_adj < isl_coord_num/2; ++isl_adj) {
            const int neighbor_i = fIslands[i].sites[j].neighbor_ind[isl_adj];
            const spins::O2 s_b = fIslands[i].sites[neighbor_i].spin;
            const spins::O2Coupling c = fIslands[i].sites[j].neighbor_coupling[isl_adj];
            total_energy -= spins::O2Product(s_a, c, s_b);
            total_energy -= spins::O2Product(s_a, fIslands[i].sites[j].local_field);
         }
      }
   }

   // TODO total_cos and total_sin
   XYChainMPI::Measures ret;
   ret.energy = total_energy;
   ret.total_spin = total_spin;
//   ret.spin_cos = total_cos;
//   ret.spin_sin = total_sin;
   ret.m_isl_avg = m_isl_avg;
   ret.m2_isl_avg = m2_isl_avg;
   return ret;
}

XYChainMPI::Measures XYChainMPI::GetMeasuresExhaustive() const {
   double total_energy = 0;
   spins::O2 total_spin = spins::O2_zero;
   vector<double> total_cos(fArrayDims.size());
   vector<double> total_sin(fArrayDims.size());
   double m_isl_avg = 0;
   double m2_isl_avg = 0;

   int arr_coord_num = fArrayCoord;
   int isl_coord_num = fIslandCoord;

   // Accumulate island total spin
   for (int i = 0; i < fNarr; ++i) {
      spins::O2 isl_total_spin = spins::O2_zero;
      for (int j = 0; j < fIslVol; ++j) {
         const spins::O2 s = fIslands[i].sites[j].spin;
         spins::O2SumAccum(s, &isl_total_spin);
      }
      fIslands[i].measures.total_spin = isl_total_spin;
      double m2 = spins::O2Product(isl_total_spin, isl_total_spin);
      spins::O2SumAccum(isl_total_spin, &total_spin);
      m_isl_avg += sqrt(m2) / fNarr;
      m2_isl_avg += m2 / fNarr;
   }

   // Accumulate inter-island coupling energy term
   for (int i = 0; i < fNarr; ++i) {
      const spins::O2 isl_total_spin = fIslands[i].measures.total_spin;
      for (int arr_adj = 0; arr_adj < arr_coord_num/2; ++arr_adj) {
         const int neighbor_i = fIslands[i].neighbor_ind[arr_adj];
         const spins::O2 isl_adj_total_spin = fIslands[neighbor_i].measures.total_spin;
         const spins::O2Coupling isl_coupling = fIslands[i].neighbor_coupling[arr_adj];
         total_energy -= spins::O2Product(isl_total_spin, isl_coupling, isl_adj_total_spin);
         total_energy -= spins::O2Product(isl_total_spin, fIslands[i].local_field);

         const spins::O2 isl_adj_rot = spins::O2Product(
               spins::O2Coupling_rotate90, isl_adj_total_spin);
//         total_cos[arr_adj] += spins::O2Product(
//               isl_total_spin, isl_coupling, isl_adj_total_spin);
//         total_sin[arr_adj] += spins::O2Product(
//               isl_total_spin, isl_coupling, isl_adj_rot);
      }
   }

   // Accumulate intra-island coupling energy terms
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         const spins::O2 s_a = fIslands[i].sites[j].spin;
         for (int isl_adj = 0; isl_adj < isl_coord_num/2; ++isl_adj) {
            const int neighbor_i = fIslands[i].sites[j].neighbor_ind[isl_adj];
            const spins::O2 s_b = fIslands[i].sites[neighbor_i].spin;
            const spins::O2Coupling c = fIslands[i].sites[j].neighbor_coupling[isl_adj];
            total_energy -= spins::O2Product(s_a, c, s_b);
            total_energy -= spins::O2Product(s_a, fIslands[i].sites[j].local_field);
         }
      }
   }

   // TODO total_cos and total_sin
   XYChainMPI::Measures ret;
   ret.energy = total_energy;
   ret.total_spin = total_spin;
//   ret.spin_cos = total_cos;
//   ret.spin_sin = total_sin;
   ret.m_isl_avg = m_isl_avg;
   ret.m2_isl_avg = m2_isl_avg;
   return ret;
}

void XYChainMPI::Sweep(void (*sweep_func)(void* arg)) {
   if (fIsDirty) {
      fMeasures = GetMeasuresExhaustive();
      fIsDirty = false;
   }

/*
   int x_max = fNslice - fNslice%2;
   for (int x = 0; x < x_max; x += 2) {
      fJobPool.AddJob(sweep_func, (void*)&fThreadArgs[x]);
   }
   fJobPool.Wait();

   x_max = fNslice;
   for (int x = 1; x < x_max; x += 2) {
      fJobPool.AddJob(sweep_func, (void*)&fThreadArgs[x]);
   }
   fJobPool.Wait();

   if (fNslice%2 == 1) {
      fJobPool.AddJob(sweep_func, (void*)&fThreadArgs[fNslice-1]);
   }
   fJobPool.Wait();
*/
   SumMeasureChanges();
}

void XYChainMPI::SumMeasureChanges() {
   fMeasures.m_isl_avg = 0;
   fMeasures.m2_isl_avg = 0;
   for (int i = 0; i < fNslice; ++i) {
      AddMeasures(fThreadArgs[i].measures_ret, &fMeasures);
   }
}

void XYChainMPI::AddMeasures(const XYChainMPI::Measures& m1, XYChainMPI::Measures* m2) {
   m2->energy += m1.energy;
   spins::O2SumAccum(m1.total_spin, &m2->total_spin);
//   for (int i = 0; i < min(m1.spin_sin.size(), m2->spin_sin.size()); ++i) {
//      m2->spin_sin[i] += m1.spin_sin[i];
//   }
//  for (int i = 0; i < min(m1.spin_cos.size(), m2->spin_cos.size()); ++i) {
//      m2->spin_cos[i] += m1.spin_cos[i];
//   }
   m2->m_isl_avg += m1.m_isl_avg;
   m2->m2_isl_avg += m1.m2_isl_avg;
}

void XYChainMPI::ZeroMeasures(XYChainMPI::Measures* m) {
   m->energy = 0;
   m->total_spin = spins::O2_zero;
//   for (int i = 0; i < m->spin_sin.size(); ++i) {
//      m->spin_sin[i] = 0;
//      m->spin_cos[i] = 0;
//   }
   m->m_isl_avg = 0;
   m->m2_isl_avg = 0;
}

void XYChainMPI::MetroSweep() {
   Sweep(MetroSweepSlice);
}

void XYChainMPI::MetroSweepSlice(void* args) {
   XYChainMPI::ThreadArgs* targs = (XYChainMPI::ThreadArgs*)args;
   XYChainMPI* self = targs->self;
   int slice_begin = targs->slice_begin;
   int slice_end = targs->slice_end;
   rng::RandomGen* rng = &targs->rng;

   ZeroMeasures(&targs->measures_ret);
   int n_arr_slice = self->fNarr / self->fArrayDims.back();
   int arr_i_begin = n_arr_slice * slice_begin;
   int arr_i_end = n_arr_slice * slice_end;

   double total_de = 0;
   spins::O2 total_ds = spins::O2_zero;
   vector<double> total_dcos;
   vector<double> total_dsin;
   double total_dm_isl = 0;
   double total_dm2_isl = 0;

   int arr_coord_num = self->fArrayCoord;
   int isl_coord_num = self->fIslandCoord;

   // Accumulate inter-island coupling energy term
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      spins::O2 arr_mf = spins::O2_zero;
      for (int arr_adj = 0; arr_adj < arr_coord_num; ++arr_adj) {
         const int neighbor_i = self->fIslands[i].neighbor_ind[arr_adj];
         const spins::O2 isl_adj_total_spin = self->fIslands[neighbor_i].measures.total_spin;
         const spins::O2Coupling isl_coupling = self->fIslands[i].neighbor_coupling[arr_adj];
         const spins::O2 isl_adj_f = spins::O2Product(isl_coupling, isl_adj_total_spin);
         spins::O2SumAccum(isl_adj_f, &arr_mf);
      }
      spins::O2SumAccum(self->fIslands[i].local_field, &arr_mf);

      for (int j = 0; j < self->fIslVol; ++j) {
         const spins::O2 s0 = self->fIslands[i].sites[j].spin;
         spins::O2 isl_mf = spins::O2_zero;
         for (int isl_adj = 0; isl_adj < isl_coord_num; ++isl_adj) {
            const int isl_neighbor_i = self->fIslands[i].sites[j].neighbor_ind[isl_adj];
            const spins::O2 s_adj = self->fIslands[i].sites[isl_neighbor_i].spin;
            const spins::O2Coupling c = self->fIslands[i].sites[j].neighbor_coupling[isl_adj];
            const spins::O2 adj_f = spins::O2Product(c, s_adj);
            spins::O2SumAccum(adj_f, &isl_mf);
         }
         spins::O2SumAccum(self->fIslands[i].sites[j].local_field, &isl_mf);

         double epsilon = targs->window * (rng->RandomUniform() * 2 - 1);
         spins::O2 s0_perp = spins::O2Product(spins::O2Coupling_rotate90, s0);
         spins::O2ScaleAccum(epsilon, &s0_perp);
         spins::O2 s1 = spins::O2Sum(s0, s0_perp);
         spins::O2ScaleAccum(1.0/sqrt(spins::O2Product(s1, s1)), &s1);

         double e0 = - spins::O2Product(s0, isl_mf) - spins::O2Product(s0, arr_mf);
         double e1 = - spins::O2Product(s1, isl_mf) - spins::O2Product(s1, arr_mf);
         double de = e1 - e0;
         if (de <= 0 || rng->RandomUniform() < exp(-de * self->fBeta)) {
            total_de += de;

            spins::O2 ds = spins::O2Diff(s1, s0);
            self->fIslands[i].sites[j].spin = s1;
            spins::O2SumAccum(ds, &self->fIslands[i].measures.total_spin);
            spins::O2SumAccum(ds, &total_ds);

            targs->n_accept++;
         } else {
            targs->n_reject++;
         }

         targs->window = (2.0 * targs->n_accept) / (targs->n_accept + targs->n_reject);
      }

      double m2_isl = spins::O2Product(self->fIslands[i].measures.total_spin,
                                      self->fIslands[i].measures.total_spin);
      total_dm_isl += sqrt(m2_isl) / self->fNarr;
      total_dm2_isl += m2_isl / self->fNarr;
   }

   targs->measures_ret.energy = total_de;
   targs->measures_ret.total_spin = total_ds;
   targs->measures_ret.m_isl_avg = total_dm_isl;
   targs->measures_ret.m2_isl_avg = total_dm2_isl;
}

void XYChainMPI::OverRelaxSweep() {
   Sweep(OverRelaxSweepSlice);
}

void XYChainMPI::OverRelaxSweepSlice(void* args) {
   XYChainMPI::ThreadArgs* targs = (XYChainMPI::ThreadArgs*)args;
   XYChainMPI* self = targs->self;
   int slice_begin = targs->slice_begin;
   int slice_end = targs->slice_end;
   rng::RandomGen* rng = &targs->rng;

   ZeroMeasures(&targs->measures_ret);
   int n_arr_slice = self->fNarr / self->fArrayDims.back();
   int arr_i_begin = n_arr_slice * slice_begin;
   int arr_i_end = n_arr_slice * slice_end;

   double total_de = 0;
   spins::O2 total_ds = spins::O2_zero;
   vector<double> total_dcos;
   vector<double> total_dsin;
   double total_dm_isl = 0;
   double total_dm2_isl = 0;

   int arr_coord_num = self->fArrayCoord;
   int isl_coord_num = self->fIslandCoord;

   // Accumulate inter-island coupling energy term
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      spins::O2 arr_mf = spins::O2_zero;
      for (int arr_adj = 0; arr_adj < arr_coord_num; ++arr_adj) {
         const int neighbor_i = self->fIslands[i].neighbor_ind[arr_adj];
         const spins::O2 isl_adj_total_spin = self->fIslands[neighbor_i].measures.total_spin;
         const spins::O2Coupling isl_coupling = self->fIslands[i].neighbor_coupling[arr_adj];
         const spins::O2 isl_adj_f = spins::O2Product(isl_coupling, isl_adj_total_spin);
         spins::O2SumAccum(isl_adj_f, &arr_mf);
      }
      spins::O2SumAccum(self->fIslands[i].local_field, &arr_mf);

      for (int j = 0; j < self->fIslVol; ++j) {
         const spins::O2 s0 = self->fIslands[i].sites[j].spin;
         spins::O2 isl_mf = spins::O2_zero;
         for (int isl_adj = 0; isl_adj < isl_coord_num; ++isl_adj) {
            const int isl_neighbor_i = self->fIslands[i].sites[j].neighbor_ind[isl_adj];
            const spins::O2 s_adj = self->fIslands[i].sites[isl_neighbor_i].spin;
            const spins::O2Coupling c = self->fIslands[i].sites[j].neighbor_coupling[isl_adj];
            const spins::O2 adj_f = spins::O2Product(c, s_adj);
            spins::O2SumAccum(adj_f, &isl_mf);
         }
         spins::O2SumAccum(self->fIslands[i].sites[j].local_field, &isl_mf);

         spins::O2 mf = spins::O2Sum(isl_mf, arr_mf);
         spins::O2 s1;
         if (spins::O2Product(mf, mf) == 0) {
            s1 = RandomUnitO2(rng);
         } else {
            s1 = spins::O2Reflect(s0, mf);
         }

         spins::O2 ds = spins::O2Diff(s1, s0);
         self->fIslands[i].sites[j].spin = s1;
         spins::O2SumAccum(ds, &self->fIslands[i].measures.total_spin);
         spins::O2SumAccum(ds, &total_ds);
      }

      double m2_isl = spins::O2Product(self->fIslands[i].measures.total_spin,
                                      self->fIslands[i].measures.total_spin);
      total_dm_isl += sqrt(m2_isl) / self->fNarr;
      total_dm2_isl += m2_isl / self->fNarr;
   }

   targs->measures_ret.energy = total_de;
   targs->measures_ret.total_spin = total_ds;
   targs->measures_ret.m_isl_avg = total_dm_isl;
   targs->measures_ret.m2_isl_avg = total_dm2_isl;
}

void XYChainMPI::FlattenStaticInfo(int arr_i_begin, int arr_i_end, char* buf) {
   int wr_ind = 0;
   *(int*)(&buf[wr_ind]) = arr_i_begin;
   wr_ind += sizeof(int);
   *(int*)(&buf[wr_ind]) = arr_i_end;
   wr_ind += sizeof(int);

   const int isl_coord_num = fIslandCoord;
   const int arr_coord_num = fArrayCoord;
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      *(spins::O2*)(&buf[wr_ind]) = fIslands[i].local_field;
      wr_ind += sizeof(spins::O2);
      for (int k = 0; k < arr_coord_num; ++k) {
         *(spins::O2Coupling*)(&buf[wr_ind]) = fIslands[i].neighbor_coupling[k];
         wr_ind += sizeof(spins::O2Coupling);
      }
      for (int k = 0; k < arr_coord_num; ++k) {
         *(int*)(&buf[wr_ind]) = fIslands[i].neighbor_ind[k];
         wr_ind += sizeof(int);
      }

      for (int j = 0; j < fIslVol; ++j) {
         *(spins::O2*)(&buf[wr_ind]) = fIslands[i].sites[j].local_field;
         wr_ind += sizeof(spins::O2);
         for (int k = 0; k < isl_coord_num; ++k) {
            *(spins::O2Coupling*)(&buf[wr_ind]) = fIslands[i].sites[j].neighbor_coupling[k];
            wr_ind += sizeof(spins::O2Coupling);
         }
         for (int k = 0; k < isl_coord_num; ++k) {
            *(int*)(&buf[wr_ind]) = fIslands[i].sites[j].neighbor_ind[k];
            wr_ind += sizeof(int);
         }
      }
   }
}

void XYChainMPI::SendStaticInfo(int arr_i_begin, int arr_i_end, int dest) {
   const int n_arr_slice = arr_i_end - arr_i_begin;
   const int n_spins = n_arr_slice * fIslVol;
   const int isl_coord_num = fIslandCoord;
   const int arr_coord_num = fArrayCoord;
   const int buflen = 1 + 2 * sizeof(int)
                        + n_arr_slice * sizeof(spins::O2)
                        + n_arr_slice * arr_coord_num * sizeof(spins::O2Coupling)
                        + n_arr_slice * arr_coord_num * sizeof(int)
                        + n_spins * sizeof(spins::O2)
                        + n_spins * isl_coord_num * sizeof(spins::O2Coupling)
                        + n_spins * isl_coord_num * sizeof(int);

   char* buf = fSendBuf;
   buf[buflen-1] = 0;

   FlattenStaticInfo(arr_i_begin, arr_i_end, buf);
   MPI::COMM_WORLD.Ssend(&buflen, 1, MPI_INT, dest, SLICE_META_BUFLEN_TAG);
   MPI::COMM_WORLD.Ssend(buf, buflen, MPI_BYTE, dest, SLICE_META_TAG);
}

void XYChainMPI::UnflattenStaticInfo(const char* buf) {
   int rd_ind = 0;
   const int arr_i_begin = *(int*)(&buf[rd_ind]);
   rd_ind += sizeof(int);
   const int arr_i_end = *(int*)(&buf[rd_ind]);
   rd_ind += sizeof(int);

   const int isl_coord_num = fIslandCoord;
   const int arr_coord_num = fArrayCoord;
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      fIslands[i].local_field = *(spins::O2*)(&buf[rd_ind]);
      rd_ind += sizeof(spins::O2);
      for (int k = 0; k < arr_coord_num; ++k) {
         fIslands[i].neighbor_coupling[k] = *(spins::O2Coupling*)(&buf[rd_ind]);
         rd_ind += sizeof(spins::O2Coupling);
      }
      for (int k = 0; k < arr_coord_num; ++k) {
         fIslands[i].neighbor_ind[k] = *(int*)(&buf[rd_ind]);
         rd_ind += sizeof(int);
      }

      for (int j = 0; j < fIslVol; ++j) {
         fIslands[i].sites[j].local_field = *(spins::O2*)(&buf[rd_ind]);
         rd_ind += sizeof(spins::O2);
         for (int k = 0; k < isl_coord_num; ++k) {
            fIslands[i].sites[j].neighbor_coupling[k] = *(spins::O2Coupling*)(&buf[rd_ind]);
            rd_ind += sizeof(spins::O2Coupling);
         }
         for (int k = 0; k < isl_coord_num; ++k) {
            fIslands[i].sites[j].neighbor_ind[k] = *(int*)(&buf[rd_ind]);
            rd_ind += sizeof(int);
         }
      }
   }

   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      for (int k = 0; k < arr_coord_num; ++k) {
         int alt_k = (k + arr_coord_num/2) % arr_coord_num;
         int alt_i = fIslands[i].neighbor_ind[k];
         fIslands[alt_i].neighbor_coupling[alt_k] = spins::O2CouplingTranspose(
               fIslands[i].neighbor_coupling[k]);
      }
   }
}

void XYChainMPI::ReceiveStaticInfo(int src) {
   int buf_len = 0;
   MPI::COMM_WORLD.Recv(&buf_len, 1, MPI_INT, src, SLICE_META_BUFLEN_TAG);
   char* buf = fSendBuf;
   buf[buf_len-1] = 0;
   MPI::COMM_WORLD.Recv(buf, buf_len, MPI_BYTE, src, SLICE_META_TAG);

   UnflattenStaticInfo(buf);
}

void XYChainMPI::FlattenSlice(int arr_i_begin, int arr_i_end, char* buf) {
   int wr_ind = 0;
   *(int*)(&buf[wr_ind]) = arr_i_begin;
   wr_ind += sizeof(int);
   *(int*)(&buf[wr_ind]) = arr_i_end;
   wr_ind += sizeof(int);

   // For each island in the slice
   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      // Need to send all the spin values, and the total measures for the island
      // Can skip things like local field, neighbor couplings since those shouldn't be dynamic
      for (int j = 0; j < fIslVol; ++j) {
         *(spins::O2*)(&buf[wr_ind]) = fIslands[i].sites[j].spin;
         wr_ind += sizeof(spins::O2);
      }
      *(XYChainMPI::Measures*)(&buf[wr_ind]) = fIslands[i].measures;
      wr_ind += sizeof(XYChainMPI::Measures);
   }
}

void XYChainMPI::SendSlice(int arr_i_begin, int arr_i_end, int dest) {
   const int n_arr_slice = arr_i_end - arr_i_begin;
   const int n_spins = n_arr_slice * fIslVol;
   const int n_measures = n_arr_slice;
   const int buflen = 1 + 2 * sizeof(int)
                        + n_spins * sizeof(spins::O2)
                        + n_measures * sizeof(XYChainMPI::Measures);

   // We'll send: slice_begin, slice_end, and then the flattened spins + measures
   char* buf = fSendBuf;
   buf[buflen-1] = 0;  // just to be safe...

   FlattenSlice(arr_i_begin, arr_i_end, buf);
   // Now just send using MPI
   MPI::COMM_WORLD.Ssend(&buflen, 1, MPI_INT, dest, SLICE_BUFLEN_TAG);
   MPI::COMM_WORLD.Ssend(buf, buflen, MPI_BYTE, dest, SLICE_TAG);
}

void XYChainMPI::UnflattenSlice(const char* buf) {
   int rd_ind = 0;
   const int arr_i_begin = *(int*)(&buf[rd_ind]);
   rd_ind += sizeof(int);
   const int arr_i_end = *(int*)(&buf[rd_ind]);
   rd_ind += sizeof(int);

   for (int i = arr_i_begin; i < arr_i_end; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         fIslands[i].sites[j].spin = *(spins::O2*)(&buf[rd_ind]);
         rd_ind += sizeof(spins::O2);
      }
      fIslands[i].measures = *(XYChainMPI::Measures*)(&buf[rd_ind]);
      rd_ind += sizeof(XYChainMPI::Measures);
   }
}

void XYChainMPI::ReceiveSlice(int src) {
   int buf_len = 0;
   MPI::COMM_WORLD.Recv(&buf_len, 1, MPI_INT, src, SLICE_BUFLEN_TAG);
   char* buf = fSendBuf;
   buf[buf_len-1] = 0;
   MPI::COMM_WORLD.Recv(buf, buf_len, MPI_BYTE, src, SLICE_TAG);

   UnflattenSlice(buf);
}

void XYChainMPI::SendMeasureDeltas(const XYChainMPI::Measures* meas) {
   const int master_id = 0;
   MPI::COMM_WORLD.Ssend(meas, sizeof(XYChainMPI::Measures), MPI_BYTE, master_id, REPORT_TAG);
   MPI::COMM_WORLD.Recv(0, 0, MPI_BYTE, master_id, REPORTS_DONE_TAG);
}

void XYChainMPI::SweepMPI(void (*sweep_func)(void* arg)) {
   const int n_ids = fNMpiIds-1; // Last ID is reserved for the master
   const int my_id = MPI::COMM_WORLD.Get_rank() - 1;
   const int down_id = 1 + (my_id + 1) % n_ids;
   const int up_id = 1 + (my_id + n_ids - 1) % n_ids;

   const int n_slice = fNarr / fArrayDims.back();
   const int a_ind = 2*my_id;
   const int b_ind = 2*my_id + 1;

   if (n_ids == 1) {
      (*sweep_func)((void*)&fThreadArgs[a_ind]);
      (*sweep_func)((void*)&fThreadArgs[b_ind]);
      XYChainMPI::Measures m = fThreadArgs[a_ind].measures_ret;
      AddMeasures(fThreadArgs[b_ind].measures_ret, &m);
      SendMeasureDeltas(&m);
      return;
   }

   (*sweep_func)((void*)&fThreadArgs[a_ind]);
   // Now the measurement deltas are stored in fThreadArgs[a_ind].measures_ret

   if (my_id%2 != 0) {
      ReceiveSlice(down_id);
   }
   // We could also probably just send the last slice instead of the entire chunk range
   SendSlice(fThreadArgs[a_ind].slice_begin * n_slice,
             (fThreadArgs[a_ind].slice_begin+1) * n_slice, up_id);
   if (my_id%2 == 0) {
      ReceiveSlice(down_id);
   }

   (*sweep_func)((void*)&fThreadArgs[b_ind]);
   // Now the deltas are in fThreadArgs[b_ind].measures_ret

   if (my_id%2 != 0) {
      ReceiveSlice(up_id);
   }
   SendSlice((fThreadArgs[b_ind].slice_end-1) * n_slice,
             fThreadArgs[b_ind].slice_end * n_slice, down_id);
   if (my_id%2 == 0) {
      ReceiveSlice(up_id);
   }

   XYChainMPI::Measures m = fThreadArgs[a_ind].measures_ret;
   AddMeasures(fThreadArgs[b_ind].measures_ret, &m);
   SendMeasureDeltas(&m);
}

void XYChainMPI::ResyncMPI() {
   const int n_ids = fNMpiIds - 1; // Last ID is reserved for the master
   const int my_id = MPI::COMM_WORLD.Get_rank() - 1;
   const int down_id = 1 + (my_id + 1) % n_ids;
   const int up_id = 1 + (my_id + n_ids - 1) % n_ids;

   const int n_slice = fNarr / fArrayDims.back();
   const int a_ind = 2*my_id;
   const int b_ind = 2*my_id + 1;

   if (n_ids == 1) {
      XYChainMPI::Measures chunk_measure_contr = GetMeasuresExhaustive();
      SendMeasureDeltas(&chunk_measure_contr);
      return;
   }

   if (my_id%2 != 0) {
      ReceiveSlice(up_id);
   }
   SendSlice((fThreadArgs[b_ind].slice_end-1) * n_slice,
             fThreadArgs[b_ind].slice_end * n_slice, down_id);
   if (my_id%2 == 0) {
      ReceiveSlice(up_id);
   }

   if (my_id%2 != 0) {
      ReceiveSlice(down_id);
   }
   SendSlice(fThreadArgs[a_ind].slice_begin * n_slice,
             (fThreadArgs[a_ind].slice_begin+1) * n_slice, up_id);
   if (my_id%2 == 0) {
      ReceiveSlice(down_id);
   }

   if (my_id%2 != 0) {
      ReceiveStaticInfo(up_id);
   }
   SendStaticInfo((fThreadArgs[b_ind].slice_end-1) * n_slice,
                  fThreadArgs[b_ind].slice_end * n_slice, down_id);
   if (my_id%2 == 0) {
      ReceiveStaticInfo(up_id);
   }

   if (my_id%2 != 0) {
      ReceiveStaticInfo(down_id);
   }
   SendStaticInfo(fThreadArgs[a_ind].slice_begin * n_slice,
                  (fThreadArgs[a_ind].slice_begin+1) * n_slice, up_id);
   if (my_id%2 == 0) {
      ReceiveStaticInfo(down_id);
   }

   XYChainMPI::Measures chunk_measure_contr = GetMeasuresExhSlice(
      fThreadArgs[a_ind].slice_begin * n_slice, fThreadArgs[b_ind].slice_end * n_slice);
   SendMeasureDeltas(&chunk_measure_contr);
}

bool XYChainMPI::DoRemoteCommand() {
   const int master_id = 0;
   int cmd_num[2];
   MPI::COMM_WORLD.Recv(cmd_num, 2, MPI_INT, master_id, COMMAND_TAG);
   const int cmd = cmd_num[0];
   const int num = cmd_num[1];

   if (cmd == COMMAND_EXIT) {
      return false;
   } else if (cmd == COMMAND_METRO_SWEEP) {
      for (int i = 0; i < num; ++i) {
         MetroSweepMPI();
      }
   } else if (cmd == COMMAND_OVERRELAX_SWEEP) {
      for (int i = 0; i < num; ++i) {
         OverRelaxSweepMPI();
      }
   } else if (cmd == COMMAND_RESYNC) {
      ResyncMPI();
   } else if (cmd == COMMAND_TEST) {
      for (int i = 0; i < num; ++i) {
         fprintf(stderr, "(%i): Rubber baby buggy bumpers\n", MPI::COMM_WORLD.Get_rank());
      }
   } else {
      fprintf(stderr, "Unknown remote command ID %i\n", cmd);
      return false;
   }
   return true;
}

void XYChainMPI::CommandLoop() {
   while(true) {
      if (!DoRemoteCommand()) break;
   }
   MPI::Finalize();
}

void XYChainMPI::MasterSendCommand(int cmd, int num) {
   const int n_ids = MPI::COMM_WORLD.Get_size();
   const int my_id = MPI::COMM_WORLD.Get_rank();
   const int master_id = 0;
   if (my_id != master_id) {
      fprintf(stderr, "MasterSendCommand called on non-master ranking node %i!\n", my_id);
      return;
   }

   int cmd_num[2];
   cmd_num[0] = cmd;
   cmd_num[1] = num;
   for (int i = 1; i < n_ids; ++i) {
      MPI::COMM_WORLD.Ssend(cmd_num, 2, MPI_INT, i, COMMAND_TAG);
   }
}

XYChainMPI::Measures XYChainMPI::MasterResync() {
   XYChainMPI::Measures result;
   ZeroMeasures(&result);

   MasterSendCommand(COMMAND_RESYNC, 0);

/*
   const int n_ids = MPI::COMM_WORLD.Get_size();
   for (int i = 1; i < n_ids; ++i) {
      XYChainMPI::Measures meas_delt;
      MPI::COMM_WORLD.Recv(&meas_delt, sizeof(XYChainMPI::Measures), MPI_BYTE, i, REPORT_TAG);
      AddMeasures(meas_delt, &result);
   }

   for (int i = 1; i < n_ids; ++i) {
      MPI::COMM_WORLD.Ssend(0, 0, MPI_BYTE, i, REPORTS_DONE_TAG);
*/

   return MasterCollect(1);
}

XYChainMPI::Measures XYChainMPI::MasterCollect(int n_sweeps) {
   XYChainMPI::Measures result;
   ZeroMeasures(&result);

   MPI::Status status;
   const int n_ids = MPI::COMM_WORLD.Get_size() - 1;
   for (int sweep = 0; sweep < n_sweeps; ++sweep) {
      for (int i = 0; i < n_ids; ++i) {
         XYChainMPI::Measures meas_delt;
         MPI::COMM_WORLD.Recv(&meas_delt, sizeof(XYChainMPI::Measures), MPI_BYTE, MPI::ANY_SOURCE, REPORT_TAG, status);
         AddMeasures(meas_delt, &result);
      }
      for (int i = 0; i < n_ids; ++i) {
         int node = i+1;
         MPI::COMM_WORLD.Ssend(0, 0, MPI_BYTE, node, REPORTS_DONE_TAG);
      }
   }

   return result;
}

XYChainMPI::Measures XYChainMPI::MasterMetroSweeps(int n_sweeps) {
   XYChainMPI::Measures result;
   ZeroMeasures(&result);

   MasterSendCommand(COMMAND_METRO_SWEEP, n_sweeps);

   return MasterCollect(n_sweeps);
}

XYChainMPI::Measures XYChainMPI::MasterOverRelaxSweeps(int n_sweeps) {
   XYChainMPI::Measures result;
   ZeroMeasures(&result);

   MasterSendCommand(COMMAND_OVERRELAX_SWEEP, n_sweeps);

   return MasterCollect(n_sweeps);
}

void XYChainMPI::MasterShutdownSlaves() {
   const int n_ids = MPI::COMM_WORLD.Get_size();
   const int my_id = MPI::COMM_WORLD.Get_rank();
   const int master_id = 0;
   if (my_id != master_id) {
      fprintf(stderr, "MasterSendCommand called on non-master ranking node %i!\n", my_id);
      return;
   }

   int cmd_num[2];
   cmd_num[0] = COMMAND_EXIT;
   cmd_num[1] = 0;
   for (int i = 1; i < n_ids; ++i) {
      MPI::COMM_WORLD.Ssend(cmd_num, 2, MPI_INT, i, COMMAND_TAG);
   }
   MPI::Finalize();
}

void XYChainMPI::MasterPrintTest() {
   MasterSendCommand(COMMAND_TEST, 3);
}

void XYChainMPI::MetroSweepMPI() {
   SweepMPI(MetroSweepSlice);
}

void XYChainMPI::OverRelaxSweepMPI() {
   SweepMPI(OverRelaxSweepSlice);
}

#ifdef VISUALS
void XYChainMPI::GLRender(double scale) {
   GLRender(0, 0, scale);
}

void XYChainMPI::GLRender(double x, double y, double scale) {
   double p_3 = atan2(0, -1) / 3;
   double p2 = atan2(0, -1) * 2;

   if (fArrayDims.size() == 1) {
      for (int i = 0; i < fArrayDims[0]; ++i) {
         double x0 = x + (i + 0.5) * scale;
         double y0 = y + 0.5 * scale;
         spins::O2 v = fIslands[i].measures.total_spin;
         spins::O2ScaleAccum(scale/2 / fIslVol, &v);

         double s = atan2(v.y, v.x);
         if (s < 0) s += p2;
         double h = s / p_3;
         double r = (h -(double)((int)h));
         if (h < 1)
            glColor3f(1, r, 0);
         else if (h < 2)
            glColor3f(1 - r, 1, 0);
         else if (h < 3)
            glColor3f(0, 1, r);
         else if (h < 4)
            glColor3f(0, 1 - r, 1);
         else if (h < 5)
            glColor3f(r, 0, 1);
         else if (h < 6)
            glColor3f(1, 0, 1 - r);

/*
         glBegin(GL_TRIANGLES);
         glVertex2f(x + i*scale, y);
         glVertex2f(x + (i+1)*scale, y);
         glVertex2f(x + (i+1)*scale, y + scale);
         glVertex2f(x + (i+1)*scale, y + scale);
         glVertex2f(x + i*scale, y + scale);
         glVertex2f(x + i*scale, y);
         glEnd();
*/
//         glColor3f(1,1,1);

         glBegin(GL_LINES);
         glVertex2f(x0 - v.x, y0 - v.y);
         glVertex2f(x0 + v.x, y0 + v.y);
         glEnd();

         glBegin(GL_TRIANGLES);
         glVertex2f(x0 + v.x, y0 + v.y);
         glVertex2f(x0 - v.y/2, y0 + v.x/2);
         glVertex2f(x0 + v.y/2, y0 - v.x/2);
         glEnd();
      }
   } else if (fArrayDims.size() >= 2) {
      for (int j = 0; j < fArrayDims[1]; ++j) {
         for (int i = 0; i < fArrayDims[0]; ++i) {
            double x0 = x + (i + 0.5) * scale;
            double y0 = y + (j + 0.5) * scale;
            int ind = j * fArrayDims[0] + i;
            spins::O2 v = fIslands[ind].measures.total_spin;
            spins::O2ScaleAccum(scale/2 / fIslVol, &v);

            double s = atan2(v.y, v.x);
            if (s < 0) s += p2;
            double h = s / p_3;
            double r = (h -(double)((int)h));
            if (h < 1)
               glColor3f(1, r, 0);
            else if (h < 2)
               glColor3f(1 - r, 1, 0);
            else if (h < 3)
               glColor3f(0, 1, r);
            else if (h < 4)
               glColor3f(0, 1 - r, 1);
            else if (h < 5)
               glColor3f(r, 0, 1);
            else if (h < 6)
               glColor3f(1, 0, 1 - r);

/*
            glBegin(GL_TRIANGLES);
            glVertex2f(x + i*scale, y + j*scale);
            glVertex2f(x + (i+1)*scale, y + j*scale);
            glVertex2f(x + (i+1)*scale, y + (j+1)*scale);
            glVertex2f(x + (i+1)*scale, y + (j+1)*scale);
            glVertex2f(x + i*scale, y + (j+1)*scale);
            glVertex2f(x + i*scale, y + j*scale);
            glEnd();
*/
//            glColor3f(1,1,1);

            glBegin(GL_LINES);
            glVertex2f(x0 - v.x, y0 - v.y);
            glVertex2f(x0 + v.x, y0 + v.y);
            glEnd();

            glBegin(GL_TRIANGLES);
            glVertex2f(x0 + v.x, y0 + v.y);
            glVertex2f(x0 - v.y/2, y0 + v.x/2);
            glVertex2f(x0 + v.y/2, y0 - v.x/2);
            glEnd();
         }
      }
   }
}
#endif  // VISUALS            
