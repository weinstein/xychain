#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"



int AnalyzeXYChain(char * datafile, Int_t eq, Int_t skip, Int_t n_step, Int_t n_sweep)
{
   float* T = new float[n_step];
   float* beta = new float[n_step];
   float* e_avg = new float[n_step];
   float* c_avg = new float[n_step];
   float* x_avg = new float[n_step];
   float* x_i_avg = new float[n_step];
   float* stiff_avg = new float[n_step];
   float* binder = new float[n_step];
   float* dbinder = new float[n_step];

   FILE *fp = fopen(datafile, "r");
   int step,sweep;
   float temp,e,m,m2,m_i,m2_i,cos_term,sin2_term;
   char line[512];

   for (int j = 0; j < n_step; ++j) {
      float e_sum,m_sum,m2_sum,m3_sum,m_i_sum,m2_i_sum,cos_term_avg,sin2_term_avg,m4_sum;
      e_sum = m_sum = m2_sum = m3_sum = m4_sum = m_i_sum = m2_i_sum 
            = cos_term_avg = sin2_term_avg = 0;
      float e2_sum = 0;
      int n_meas = 0;
      for (int i = 0; i < n_sweep; ++i) {
        fgets(&line,512,fp);
      	sscanf(&line[0],"%i %i %f %f %f %f %f %f %f %f",
               &step,&sweep,&temp,&e,&m,&m2,&m_i,&m2_i,&cos_term,&sin2_term);
        if (sweep < eq) continue;
      	if (sweep%(skip+1) > 0) continue;
      	e_sum += e;
        e2_sum += e * e;
      	m_sum += m;
      	m2_sum += m2;
        m3_sum += m*m*m;
        m4_sum += m*m*m*m;
      	m_i_sum += m_i;
      	m2_i_sum += m2_i;
        cos_term_avg += cos_term;
        sin2_term_avg += sin2_term;
      	n_meas++;
      }
      T[step] = temp;
      beta[step] = 1.0/temp;
      e_avg[step] = e_sum/n_meas;
      float e2_avg = e2_sum/n_meas;
      c_avg[step] = e2_sum - e_sum*e_sum;
      float mavg = m_sum/n_meas;
      float m2avg = m2_sum/n_meas;
      float m3avg = m3_sum/n_meas;
      float m4avg = m4_sum/n_meas;
      binder[step] = (m4avg - 4*mavg*m3avg + 6*mavg*mavg*m2avg - 3*mavg*mavg*mavg*mavg)
            / ((m2avg - mavg*mavg)*(m2avg - mavg*mavg));
      if (step >= 4) {
         dbinder[step-2] = (-binder[step] + 8*binder[step-1] - 8*binder[step-3] + binder[step-4])
               / (3 * (beta[step] - beta[step-4]));
      } else {
         if (step == 1) dbinder[0] = (binder[1] - binder[0]) / (beta[1] - beta[0]);
         else if (step == 2) {
            dbinder[step-1] = (binder[step] - binder[step-2]) / (beta[step] - beta[step-2]);
         }
         else {
            dbinder[step] = (binder[step] - binder[step-1]) / (beta[step] - beta[step-1]);
         }
      }
      if (step == n_step-1 || step == n_step-2) {
         dbinder[step] = (binder[step] - binder[step-1]) / (beta[step] - beta[step-1]);
      }

      x_avg[step] = m2avg - (mavg*mavg);
      float miavg = m_i_sum/n_meas;
      float m2iavg = m2_i_sum/n_meas;
      x_i_avg[step] = m2iavg - (miavg*miavg);
      float ct = cos_term_avg/n_meas;
      float st2 = sin2_term_avg/n_meas;
      stiff_avg[step] = ct - 1.0/temp * st2;
   }

   for (int i = 0; i < n_step; ++i) {
      printf("%i %f %f %f %f %f %f %f %f %f\n",i,T[i],beta[i],e_avg[i],x_avg[i],x_i_avg[i],
            stiff_avg[i],binder[i],dbinder[i],c_avg[i]);
   }
}
