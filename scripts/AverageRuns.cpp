#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char ** argv)
{
   std::vector<FILE*> fps;
   for (int i = 1; i < argc; ++i) {
      fps.push_back(fopen(argv[i], "r"));
      if (fps.back() == 0) {
         printf("Error opening file %s\n", argv[i]);
         return 1;
      }
   }

   const int n_files = fps.size();
   char line[512];
   while (true) {
      float e_avg, e2_avg, m_avg, m2_avg, m_i_avg, m2_i_avg;
      e_avg = e2_avg = m_avg = m2_avg = m_i_avg = m2_i_avg = 0;

      int step, sweep;
      step = sweep = 0;
      float temp = 0;

      for (int i = 0; i < n_files; ++i) {
         float e, m, m2, m_i, m2_i;
         e = m = m2 = m_i = m2_i = 0;
         fgets(&line[0], 512, fps[i]);
         if (feof(fps[i])) {
            return 0;
         }
         sscanf(&line[0], "%i %i %f %f %f %f %f %f", &step, &sweep, &temp, &e, &m, &m2, &m_i, &m2_i);
         e_avg += e;
         e2_avg += e*e;
         m_avg += m;
         m2_avg += m2;
         m_i_avg += m_i;
         m2_i_avg += m2_i;
      }
      e_avg /= n_files;
      e2_avg /= n_files;
      m_avg /= n_files;
      m2_avg /= n_files;
      m_i_avg /= n_files;
      m2_i_avg /= n_files;
      printf("%i %i %f %f %f %f %f %f %f %f %f\n", step, sweep, temp, e_avg, m_avg, m2_avg, m_i_avg, m2_i_avg, sqrt(e2_avg - e_avg * e_avg), sqrt(m2_avg - m_avg * m_avg), sqrt(m2_i_avg - m_i_avg * m_i_avg));
   }
}
