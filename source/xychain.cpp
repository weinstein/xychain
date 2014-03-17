#include "xychain.h"
#include <math.h>
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
using std::endl;

#include "notification.h"

#ifdef VISUALS
 #ifdef LINUX
  #include <GL/gl.h>
 #endif  // LINUX

 #ifdef OSX
  #include <OpenGL/gl.h>
 #endif  // OSX
#endif  // VISUALS

XYChain::XYChain(vector<int> array_dimensions,
                 vector<int> island_dimensions,
                 int n_threads)
      : fNthreads(n_threads),
        fJobPool(fNthreads),
        fArrayDims(array_dimensions),
        fIslandDims(island_dimensions),
        fBeta(1),
        fNaccept(0),
        fNreject(0),
        fWindow(1) {
   InitIslands();
   InitThreading();
   fIsDirty = true;
}

void XYChain::InitIslands() {
   fNarr = 1;
   for (int i = 0; i < fArrayDims.size(); ++i) {
      fNarr *= fArrayDims[i];
   }

   fIslVol = 1;
   for (int i = 0; i < fIslandDims.size(); ++i) {
      fIslVol *= fIslandDims[i];
   }

   fTotalVol = fNarr * fIslVol;

   fIslands = new XYChain::Island[fNarr];
   int arr_coord_num = 2*fArrayDims.size();
   int isl_coord_num = 2*fIslandDims.size();
   for (int i = 0; i < fNarr; ++i) {
      fIslands[i].sites = new XYChain::SpinSite[fIslVol];
      for (int j = 0; j < fIslVol; ++j) {
         fIslands[i].sites[j].neighbor_coupling = new spins::O2Coupling[isl_coord_num];
         fIslands[i].sites[j].neighbor_ind = MakeNeighborIndices(fIslandDims, j);
         fIslands[i].sites[j].local_field = {0, 0};
      }
      fIslands[i].neighbor_coupling = new spins::O2Coupling[arr_coord_num];
      fIslands[i].neighbor_ind = MakeNeighborIndices(fArrayDims, i);
      fIslands[i].local_field = {0, 0};
   }

   fIsDirty = true;
}

void XYChain::InitThreading() {
   fNslice = fNthreads * 2;
   if (fArrayDims.back() < fNslice) {
      fNslice = fArrayDims.back();
   }
   fThreadArgs = new XYChain::ThreadArgs[fNslice];
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

void XYChain::ClearThreading() {
   delete [] fThreadArgs;
}

void XYChain::ClearIslands() {
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

XYChain::~XYChain() {
   ClearIslands();
   ClearThreading();
}

int* XYChain::MakeNeighborIndices(const std::vector<int>& dims, int index) {
   int coord_num = 2*dims.size();
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
      int neighbor_index = 0;
      int pow = 1;
      // Convert coordinates back into flat index
      for (int dim_j = 0; dim_j < dims.size(); ++dim_j) {
         neighbor_index += pow * index_coords[dim_j];
         pow *= dims[dim_j];
      }
      neighbor_i[dim_i] = neighbor_index;
      index_coords[dim_i] = (index_coords[dim_i]+dims[dim_i]-1) % dims[dim_i];
   }
   // Do it again, for negative neighbors
   for (int dim_i = 0; dim_i < dims.size(); ++dim_i) {
      index_coords[dim_i] = (index_coords[dim_i]+dims[dim_i]-1) % dims[dim_i];
      int neighbor_index = 0;
      int pow = 1;
      // Convert coordinates back into flat index
      for (int dim_j = 0; dim_j < dims.size(); ++dim_j) {
         neighbor_index += pow * index_coords[dim_j];
         pow *= dims[dim_j];
      }
      neighbor_i[dim_i+dims.size()] = neighbor_index;
      index_coords[dim_i] = (index_coords[dim_i]+1) % dims[dim_i];
   }
   // We're done!
   return neighbor_i;
}

spins::O2 XYChain::RandomUnitO2(rng::RandomGen* rng) {
   double pi = acos(0)*2;
   double angle = rng->RandomUniform()*2*pi;
   double x = cos(angle);
   double y = sin(angle);
   return {x,y};
}

void XYChain::Init(XYChain::Config initial_config) {
   rng::RandomGen rng;
   int arr_coord_num = 2*fArrayDims.size();
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
         int isl_coord_num = 2*fIslandDims.size();
         for (int coupl_i = 0; coupl_i < isl_coord_num; ++coupl_i) {
            fIslands[i].sites[site_i].neighbor_coupling[coupl_i] =
               spins::O2Coupling_I;
         }
      }
      fIslands[i].measures.total_spin = total_spin;
   }

   fIsDirty = true;
}

void XYChain::ScaleArrayCouplingAmplitude(double r) {
   int arr_coord_num = 2*fArrayDims.size();
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < arr_coord_num; ++j) {
         spins::O2ScaleCouplingAccum(r, &fIslands[i].neighbor_coupling[j]);
      }
   }

   fIsDirty = true;
}

void XYChain::ScaleIslandCouplingAmplitude(double r) {
   int isl_coord_num = 2*fIslandDims.size();
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

void XYChain::InitArrayCoupling(double r, XYChain::Disorder disorder) {
   rng::RandomGen rng;
   int arr_coord_num = 2*fArrayDims.size();
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

void XYChain::InitArrayCoupling(double r, XYChain::Disorder disorder, double param) {
   rng::RandomGen rng;
   int arr_coord_num = 2*fArrayDims.size();
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

void XYChain::InitIslandField(double r, XYChain::Config type) {
   rng::RandomGen rng;
   int isl_coord_num = 2*fIslandDims.size();
   for (int i = 0; i < fNarr; ++i) {
      for (int j = 0; j < fIslVol; ++j) {
         spins::O2 f = spins::O2_i;
         if (type == XYChain::CONFIG_RANDOM) {
            f = RandomUnitO2(&rng);
         }
         spins::O2ScaleAccum(r, &f);
         fIslands[i].sites[j].local_field = f;
      }
   }

   fIsDirty = true;
}

void XYChain::InitArrayField(double r, XYChain::Config type) {
   rng::RandomGen rng;
   int arr_coord_num = 2*fArrayDims.size();
   for (int i = 0; i < fNarr; ++i) {
      spins::O2 f = spins::O2_i;
      if (type == XYChain::CONFIG_RANDOM) {
         f = RandomUnitO2(&rng);
      }
      spins::O2ScaleAccum(r, &f);
      fIslands[i].local_field = f;
   }

   fIsDirty = true;
}

void XYChain::InitIslandCoupling(double r, XYChain::Disorder disorder) {
   rng::RandomGen rng;
   int isl_coord_num = 2*fIslandDims.size();
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

void XYChain::InitIslandCoupling(double r, XYChain::Disorder disorder, double param) {
   rng::RandomGen rng;
   int isl_coord_num = 2*fIslandDims.size();
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

spins::O2Coupling XYChain::RandomCoupling(rng::RandomGen* rng, double r, Disorder type) {
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

spins::O2Coupling XYChain::RandomCoupling(rng::RandomGen* rng, double r, Disorder type,
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

XYChain::Measures XYChain::GetMeasuresExhaustive() const {
   double total_energy = 0;
   spins::O2 total_spin = spins::O2_zero;
   vector<double> total_cos(fArrayDims.size());
   vector<double> total_sin(fArrayDims.size());
   double m_isl_avg = 0;
   double m2_isl_avg = 0;

   int arr_coord_num = 2*fArrayDims.size();
   int isl_coord_num = 2*fIslandDims.size();

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
         total_cos[arr_adj] += spins::O2Product(
               isl_total_spin, isl_coupling, isl_adj_total_spin);
         total_sin[arr_adj] += spins::O2Product(
               isl_total_spin, isl_coupling, isl_adj_rot);
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
   XYChain::Measures ret;
   ret.energy = total_energy;
   ret.total_spin = total_spin;
   ret.spin_cos = total_cos;
   ret.spin_sin = total_sin;
   ret.m_isl_avg = m_isl_avg;
   ret.m2_isl_avg = m2_isl_avg;
   return ret;
}

void XYChain::Sweep(void (*sweep_func)(void* arg)) {
   if (fIsDirty) {
      fMeasures = GetMeasuresExhaustive();
      fIsDirty = false;
   }

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

   SumMeasureChanges();
}

void XYChain::SumMeasureChanges() {
   fMeasures.m_isl_avg = 0;
   fMeasures.m2_isl_avg = 0;
   for (int i = 0; i < fNslice; ++i) {
      AddMeasures(fThreadArgs[i].measures_ret, &fMeasures);
   }
}

void XYChain::AddMeasures(const XYChain::Measures& m1, XYChain::Measures* m2) {
   m2->energy += m1.energy;
   spins::O2SumAccum(m1.total_spin, &m2->total_spin);
   for (int i = 0; i < m1.spin_sin.size(); ++i) {
      m2->spin_sin[i] += m1.spin_sin[i];
      m2->spin_cos[i] += m1.spin_cos[i];
   }
   m2->m_isl_avg += m1.m_isl_avg;
   m2->m2_isl_avg += m1.m2_isl_avg;
}

void XYChain::ZeroMeasures(XYChain::Measures* m) {
   m->energy = 0;
   m->total_spin = spins::O2_zero;
   for (int i = 0; i < m->spin_sin.size(); ++i) {
      m->spin_sin[i] = 0;
      m->spin_cos[i] = 0;
   }
   m->m_isl_avg = 0;
   m->m2_isl_avg = 0;
}

void XYChain::MetroSweep() {
   Sweep(MetroSweepSlice);
}

void XYChain::MetroSweepSlice(void* args) {
   XYChain::ThreadArgs* targs = (XYChain::ThreadArgs*)args;
   XYChain* self = targs->self;
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

   int arr_coord_num = 2*self->fArrayDims.size();
   int isl_coord_num = 2*self->fIslandDims.size();

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

void XYChain::OverRelaxSweep() {
   Sweep(OverRelaxSweepSlice);
}

void XYChain::OverRelaxSweepSlice(void* args) {
   XYChain::ThreadArgs* targs = (XYChain::ThreadArgs*)args;
   XYChain* self = targs->self;
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

   int arr_coord_num = 2*self->fArrayDims.size();
   int isl_coord_num = 2*self->fIslandDims.size();

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

#ifdef VISUALS
void XYChain::GLRender(double scale) {
   GLRender(0, 0, scale);
}

void XYChain::GLRender(double x, double y, double scale) {
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
