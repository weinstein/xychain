#ifndef SPINS_H
#define SPINS_H

#include <iostream>
using std::ostream;
#include <string.h>

namespace spins {

struct O2 {
   double x;
   double y;
};

struct O2Coupling {
   double a, b;
   double c ,d;
};

const O2 O2_zero = {0,0};
const O2 O2_i = {1,0};
const O2 O2_j = {0,1};

const O2Coupling O2Coupling_zero = {0,0,0,0};
const O2Coupling O2Coupling_I = {1,0,0,1};
const O2Coupling O2Coupling_negI = {-1,0,0,-1};
const O2Coupling O2Coupling_rotate90 = {0, 1, -1, 0};

inline ostream& operator<<(ostream& stream, const O2& s) {
   return stream << "(" << s.x << "," << s.y << ")";
}

inline ostream& operator<<(ostream& stream, const O2Coupling& c) {
   return stream << "(" << c.a << "," << c.b << ";" << c.c << "," << c.d << ")";
}

inline double O2Product(const O2& s1, const O2& s2) {
   return s1.x*s2.x + s1.y*s2.y;
}

inline void O2SumAccum(const O2& s1, O2* s2) {
   s2->x += s1.x;
   s2->y += s1.y;
}

inline O2 O2Sum(const O2& s1, const O2& s2) {
   return {s1.x + s2.x, s1.y + s2.y};
}

inline void O2ProductAccum(const O2Coupling& c, O2* s) {
   double tx = s->x;
   double ty = s->y;
   s->x = c.a * tx + c.b * ty;
   s->y = c.c * tx + c.d * ty;
}

inline O2 O2Product(const O2Coupling& c, const O2& s) {
   return {c.a * s.x + c.b * s.y, c.c * s.x + c.d * s.y};
}

inline double O2Product(const O2& s1,
                        const O2Coupling& c,
                        const O2& s2) {
   return (s1.x*(c.a*s2.x + c.b*s2.y) + s1.y*(c.c*s2.x + c.d*s2.y));
}

inline O2 O2Scale(const double r, const O2& s) {
   return {r*s.x, r*s.y};
}

inline void O2ScaleAccum(const double r, O2* s) {
   s->x *= r;
   s->y *= r;
}

inline O2Coupling O2ScaleCoupling(const double r, const O2Coupling& c) {
   return {r*c.a, r*c.b, r*c.c, r*c.d};
}

inline void O2ScaleCouplingAccum(const double r, O2Coupling* c) {
   c->a *= r;
   c->b *= r;
   c->c *= r;
   c->d *= r;
}

inline O2 O2Diff(const O2& s1, const O2& s2) {
   return {s1.x - s2.x, s1.y - s2.y};
}

inline void O2DiffAccum(const O2& s1, O2* s2) {
   s2->x = s1.x - s2->x;
   s2->y = s1.y - s2->y;
}

inline O2 O2Reflect(const O2& s, const O2& axis) {
   double r = 2 * (s.x * axis.x + s.y * axis.y) / (axis.x*axis.x + axis.y*axis.y);
   return {-s.x + r*axis.x, -s.y + r*axis.y};
}

inline O2Coupling O2CouplingTranspose(const O2Coupling& c) {
   return {c.a, c.c, c.b, c.d};
}

inline bool operator==(const O2& lhs, const O2& rhs) {
   return (memcmp(&lhs, &rhs, sizeof(O2)) == 0);
}

inline bool operator==(const O2Coupling& lhs, const O2Coupling& rhs) {
   return (memcmp(&lhs, &rhs, sizeof(O2Coupling)) == 0);
}

}  // namespace spins

#endif
