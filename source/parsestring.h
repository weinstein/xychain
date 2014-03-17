#ifndef PARSESTRING_H
#define PARSESTRING_H

#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>

inline bool HasPrefix(const std::string& s, const std::string& pref) {
   return s.compare(0, pref.size(), pref) == 0;
}

inline std::vector<std::string> split(const std::string& str, char delim) {
   std::vector<std::string> elems;
   std::stringstream ss(str);
   std::string item;
   while (std::getline(ss, item, delim)) {
      elems.push_back(item);
   }
   return elems;
}

inline std::vector<int> ParseIntList(const std::string& str) {
   std::vector<int> ret;
   std::vector<std::string> elems = split(str, ',');
   for (const std::string& item : elems) {
      ret.push_back(atoi(item.c_str()));
   }
   return ret;
}

inline std::vector<double> ParseDoubleList(const std::string& str) {
   std::vector<double> ret;
   std::vector<std::string> elems = split(str, ',');
   for (const std::string& item : elems) {
      ret.push_back(atof(item.c_str()));
   }
   return ret;
}

// If the string has a parameter, return true and store the parameter value into *result.
// For example, ParseDoubleParam on "xyz[1.23]" returns true with *result = 1.23
// ParseDoubleParam on "xyz" returns false
inline bool ParseDoubleParam(const std::string& str, double* result) {
   std::size_t pos0 = str.find("[");
   std::size_t pos1 = str.find("]");
   if (pos0 == std::string::npos || pos1 == std::string::npos) return false;
   std::string sub = str.substr(pos0+1, pos1-pos0-1);
   *result = atof(sub.c_str());
   return true;
}

#endif
