#ifndef TESTING_H
#define TESTING_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
using std::cout;
using std::endl;
using std::string;

#define DEFINE_TEST(name) \
void  __do_##name(); \
struct name { \
   name() { \
      cout << string(40, '=') << endl; \
      cout << __func__ << "..." << endl;\
      __do_##name(); \
      cout << __func__ << " done" << endl; \
      cout << string(40, '=') << endl << endl; \
   } \
}; \
name __testing_##name; \
void __do_##name()

#define __ASSERT(lhs, rhs, cmp) \
   if ((lhs) cmp (rhs)); \
   else { \
      cout << __FILE__ << ":" << __LINE__ << endl;\
      cout << #lhs << " (" << (lhs) << ")"\
      << " " << #cmp << " " \
      << #rhs << " (" << (rhs) << ")" << endl \
      << "Assertion failed." << endl; \
      exit(1); \
   }

#define ASSERT_EQ(lhs, rhs) __ASSERT(lhs, rhs, ==)
#define ASSERT_NEQ(lhs, rhs) __ASSERT(lhs, rhs, !=)
#define ASSERT_LT(lhs, rhs) __ASSERT(lhs, rhs, <)
#define ASSERT_GT(lhs, rhs) __ASSERT(lhs, rhs, >)
#define ASSERT_LTEQ(lhs, rhs) __ASSERT(lhs, rhs, <=)
#define ASSERT_GTEQ(lhs, rhs) __ASSERT(lhs, rhs, >=)

#endif
