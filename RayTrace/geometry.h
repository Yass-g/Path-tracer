#include <math.h>   
#include <stdio.h>  
#include <cstdlib>
#include <stdlib.h>
#include <random>

#define M_PI 3.1415926
const double eps = 1e-6;
using namespace std;

double bound(double x) {
  return x < 0 ? 0 : x > 1 ? 1 : x;
}
int toInt(double x) {
  return int(pow(bound(x), 1 / 2.2) * 255 + .5);
}
double erand48(int X) {
  return rand() / double(RAND_MAX);

}

enum ref_typ {
  DIFF,
  SPEC,
  REFR
}; // Material type : diffuse, specular, refractive
// TO DO : add anisotropic and texture handling	
struct Vect {
  double x, y, z;
  Vect(double mx = 0, double my = 0, double mz = 0) {
    x = mx;
    y = my;
    z = mz;
  }

  Vect operator + (const Vect & b) const {
    return Vect(x + b.x, y + b.y, z + b.z);
  }
  Vect operator - (const Vect & b) const {
    return Vect(x - b.x, y - b.y, z - b.z);
  }
  Vect operator * (double b) const {
    return Vect(x * b, y * b, z * b);
  }

  Vect mult(const Vect & b) const {
    return Vect(x * b.x, y * b.y, z * b.z);
  }
  Vect & norm() {
    return *this = * this * (1 / sqrt(x * x + y * y + z * z));
  }
  double dot(const Vect & b) const {
    return x * b.x + y * b.y + z * b.z;
  }
  Vect operator % (Vect & b) {
    return Vect(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
  } //cross product
};
struct Ray { //Ray  structure
  Vect o, d; // origin and direction
  Vect invdir;
  int sgn[3];
  Ray(Vect ori, Vect dir): o(ori), d(dir) {
    invdir.x = 1 / d.x;
    invdir.y = 1 / d.y;
    invdir.z = 1 / d.z;
    sgn[0] = (invdir.x < 0);
    sgn[1] = (invdir.y < 0);
    sgn[2] = (invdir.z < 0);
  }
};

class Object {
  public:
    Vect e, c; // emission, color
  ref_typ refl; // Material
  virtual double intersect(const Ray & r) const = 0;
  virtual Vect normal(Vect & x) const = 0;
};
class Sphere: public Object {
  public: double rad; // radius
  Vect p; // position
  Sphere(double radi, Vect pos, Vect em, Vect col, ref_typ reflec) {
    rad = radi;
    p = pos;
    e = em;
    c = col;
    refl = reflec;
  }
  Vect normal(Vect & x) const {
    return (x - p).norm();
  }
  double intersect(const Ray & r) const { // distance of intersection (0 if not)
    Vect op = p - r.o; // Solve the line/sphere intersection equation
    double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
    if (det < 0) return 0;
    else det = sqrt(det);
    return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
  }
};

class Cube: public Object {

  public: Vect bounds[2];
  Cube(double d, Vect b1, Vect em, Vect col, ref_typ reflec) {
    bounds[0] = Vect(b1.x - d, b1.y - d, b1.z - d);
    bounds[1] = Vect(b1.x + d, b1.y + d, b1.z + d);
    e = em;
    c = col;
    refl = reflec;

  }
  Vect normal(Vect & x) const {
    Vect C = Vect((bounds[0].x + bounds[1].x) / 2, (bounds[0].y + bounds[1].y) / 2, (bounds[0].z + bounds[1].z) / 2);
    Vect p = x - C;
    Vect d = Vect(abs(bounds[0].x + bounds[1].x) / 2, abs(bounds[0].y + bounds[1].y) / 2, abs(bounds[0].z + bounds[1].z) / 2);
    Vect n = Vect(int(p.x / d.x), int(p.y / d.y), int(p.z / d.z));
    return n.norm();
  }
  double intersect(const Ray & r) const {
    double t;
    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[r.sgn[0]].x - r.o.x) * r.invdir.x;
    tmax = (bounds[1 - r.sgn[0]].x - r.o.x) * r.invdir.x;
    tymin = (bounds[r.sgn[1]].y - r.o.y) * r.invdir.y;
    tymax = (bounds[1 - r.sgn[1]].y - r.o.y) * r.invdir.y;

    if ((tmin > tymax) || (tymin > tmax))
      return 0;

    if (tymin > tmin)
      tmin = tymin;
    if (tymax < tmax)
      tmax = tymax;

    tzmin = (bounds[r.sgn[2]].z - r.o.z) * r.invdir.z;
    tzmax = (bounds[1 - r.sgn[2]].z - r.o.z) * r.invdir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
      return 0;

    if (tzmin > tmin)
      tmin = tzmin;
    if (tzmax < tmax)
      tmax = tzmax;

    t = tmin;

    if (t < 0) {
      t = tmax;
      if (t < 0) return 0;
    }

    return t;
  }

};

class Plane: public Object {
  public: Vect n;
  double d;
  Plane(double d_, Vect n_, Vect e_, Vect c_, ref_typ reflec) {
    d = d_;
    n = n_;
    e = e_;
    c = c_;
    refl = reflec;
  }
  double intersect(const Ray & ray) const {
    double d0 = n.dot(ray.d);
    if (d0 != 0) {
      double t = -1 * (((n.dot(ray.o)) + d) / d0);
      return (t > eps) ? t : 0;
    } else return 0;
  }
  Vect normal(const Vect & p0) const {
    return n;
  }
};
