#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
using std::string;

int annealing2ascii(char* fname);
int hyst2ascii(char* fname);

int main(int argc, char** argv) {
   if (argc != 3) {
      std::cerr << "Invalid argc\n";
      std::cerr << "Usage:\n";
      std::cerr << "    " << argv[0] << " [annealing|hyst] inputfile\n\n";
      return -1;
   }

   // Output format may depend on which command produced the output
   if (string(argv[1]) == "annealing") {
      return annealing2ascii(argv[2]);
   } else if (string(argv[1]) == "hyst") {
      return hyst2ascii(argv[2]);
   } else {
      std::cerr << "Invalid: " << argv[1] << "\n";
      return -1;
   }
}

int annealing2ascii(char* fname) {
   FILE* fp = fopen(fname, "r");
   while (!feof(fp)) {
      int step, sweep;
      double temp, e, m, m2, m_isl, m2_isl, mx, my;
      fread(&step, sizeof(int), 1, fp);
      fread(&sweep, sizeof(int), 1, fp);
      fread(&temp, sizeof(double), 1, fp);
      fread(&e, sizeof(double), 1, fp);
      fread(&m, sizeof(double), 1, fp);
      fread(&m2, sizeof(double), 1, fp);
      fread(&m_isl, sizeof(double), 1, fp);
      fread(&m2_isl, sizeof(double), 1, fp);
      fread(&mx, sizeof(double), 1, fp);
      fread(&my, sizeof(double), 1, fp);

      printf("%i %i %.16f ", step, sweep, temp);
      printf("%.16f %.16f %.16f ", e, m, m2);
      printf("%.16f %.16f ", m_isl, m2_isl);
      printf("%.16f %.16f ", mx, my);
      printf("\n");
   }

   return 0;
}

int hyst2ascii(char* fname) {
   FILE* fp = fopen(fname, "r");
   while (!feof(fp)) {
      int cycle, step, sweep;
      double field, e, m, m2, m_isl, m2_isl, mx, my;
      fread(&cycle, sizeof(int), 1, fp);
      fread(&step, sizeof(int), 1, fp);
      fread(&sweep, sizeof(int), 1, fp);
      fread(&field, sizeof(double), 1, fp);
      fread(&e, sizeof(double), 1, fp);
      fread(&m, sizeof(double), 1, fp);
      fread(&m2, sizeof(double), 1, fp);
      fread(&m_isl, sizeof(double), 1, fp);
      fread(&m2_isl, sizeof(double), 1, fp);
      fread(&mx, sizeof(double), 1, fp);
      fread(&my, sizeof(double), 1, fp);

      printf("%i %i %i %.16f ", cycle, step, sweep, field);
      printf("%.16f %.16f %.16f ", e, m, m2);
      printf("%.16f %.16f ", m_isl, m2_isl);
      printf("%.16f %.16f ", mx, my);
      printf("\n");
   }

   return 0;
}
