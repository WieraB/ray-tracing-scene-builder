#ifndef INTERSECTABLE_H
#define INTERSECTABLE_H

#include <math.h> 
#include <stdlib.h>
#include <stdio.h> 
#include "../Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "ray.h"

#define EIGEN_NO_DEBUG


// General base class for all objects in the scene
// Emission, colour, and reflection type should be specified for all objects.
// For the specific ray, intersect function calculates the distance to the object
// and the normal vector at the intersection point. 
class intersectable {

  public:
    Eigen::Vector3d e, c;      // emission, colour
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    virtual ~intersectable() = default;

    intersectable(Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
    e(e_), c(c_), refl(refl_) {};
    virtual double intersect(const Ray &r, Eigen::Vector3d &n_out) const = 0;

};

// Additional radius and centre poistion parameters
class Sphere : public intersectable {

  public:
    double rad;       // radius
    Eigen::Vector3d p;      // position
    Sphere(double rad_, Eigen::Vector3d p_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
    intersectable(e_, c_, refl_), rad(rad_), p(p_) {}

    double intersect(const Ray &r, Eigen::Vector3d &n_out) const override { // returns distance, 0 if nohit
      Eigen::Vector3d op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
      double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
      if (det<0) return 0; else det=sqrt(det);
  
      Eigen::Vector3d x;
      t = b - det;
      x = r.o + r.d * t;
      if (t > eps) {
          n_out = (x - p).normalized();
          return t;
      }
  
      t = b + det;
      x = r.o + r.d * t;
      if (t > eps) {
          n_out = (x - p).normalized();
          return t;
      }
      return 0;
    }

};

// Additional coorinates of 3 triangle verticies
class Triangle : public intersectable {

  public:
    Eigen::Vector3d v0, v1, v2;
  
    Triangle(Eigen::Vector3d v0_, Eigen::Vector3d v1_, Eigen::Vector3d v2_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
      intersectable(e_, c_, refl_), v0(v0_), v1(v1_), v2(v2_) {}
  
    double intersect(const Ray &r, Eigen::Vector3d &n_out) const {
      const double eps = 1e-5;
      Eigen::Vector3d edge1 = v1 - v0, edge2 = v2 - v0;
      Eigen::Vector3d pvec = r.d.cross(edge2);
      double det = edge1.dot(pvec);
      if (fabs(det) < eps) return 0;
      double invDet = 1.0 / det;
      Eigen::Vector3d tvec = r.o - v0;
      double u = tvec.dot(pvec) * invDet;
      if (u < 0 || u > 1) return 0;
      Eigen::Vector3d qvec = tvec.cross(edge1);
      double v = r.d.dot(qvec) * invDet;
      if (v < 0 || u + v > 1) return 0;
      double t = edge2.dot(qvec) * invDet;
      if (t < eps) return 0;
      n_out = edge1.cross(edge2);
      n_out.normalize();
      return t;
    }
};

std::vector<Sphere> spheres;

std::vector<Triangle> triangles;

inline bool intersect(const Ray &r, double &t, int &id, bool &isTri, Eigen::Vector3d &n,
std::vector<Sphere> &spheres, std::vector<Triangle> triangles) {
  double inf=t=1e20;
  isTri=false;
  double d;
  for (int i=0; i<(int)spheres.size(); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=spheres[i].intersect(r, n_tmp)) && d<t) { t=d; id=i; isTri=false; n=n_tmp; }
  }
  for (int i=0; i<(int)triangles.size(); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=triangles[i].intersect(r, n_tmp)) && d<t) {
      t=d; id=i; isTri=true; n=n_tmp;
    }
  }
  return t<inf;
}

Eigen::Vector3d radiance(const Ray &r, int depth, unsigned short *Xi){

  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  bool isTri=false;
  Eigen::Vector3d n;

  if (!intersect(r, t, id, isTri, n, spheres, triangles)) return Eigen::Vector3d(0, 0, 0); // if miss, return black
  Eigen::Vector3d x = r.o + r.d * t;
  Eigen::Vector3d nl = n.dot(r.d) < 0 ? n : n * -1; // Adjust normal direction
  Eigen::Vector3d f, e;
  Refl_t refl;

  if (!isTri) {
    const Sphere &obj=spheres[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  } else {
    const Triangle &obj=triangles[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  }
  
  double p = f.x()>f.y() && f.x()>f.z() ? f.x() : f.y()>f.z() ? f.y() : f.z(); // max refl

  if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return e; //R.R.

  if (refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Eigen::Vector3d w=nl;

    Eigen::Vector3d u;
    if (std::fabs(w.x()) > 0.1) {
        u = Eigen::Vector3d(0, 1, 0).cross(w);
    } else {
        u = Eigen::Vector3d(1, 0, 0).cross(w);
    }
    u.normalize();  // Normalizing the vector
    
    Eigen::Vector3d v = w.cross(u);

    Eigen::Vector3d d = u * std::cos(r1) * r2s + v * std::sin(r1) * r2s + w * std::sqrt(1 - r2);
    d.normalize();

    return e + f.cwiseProduct(radiance(Ray(x,d),depth,Xi));

  } 
  else if (refl == SPEC)            // Ideal SPECULAR reflection
    return e + f.cwiseProduct(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));

  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
  
    bool into = n.dot(nl)>0;                // Ray from outside going in?

  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;

  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return e + f.cwiseProduct(radiance(reflRay,depth,Xi));
  
  Eigen::Vector3d tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + std::sqrt(cos2t))));
  tdir.normalize();

  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
  
  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);

    // Russian roulette for depth > 2
  if (depth > 2) {
        if (erand48(Xi) < P) {
            // If we continue with the ray, calculate reflected radiance
            return e + f.cwiseProduct(radiance(reflRay, depth, Xi) * RP);
        } else {
            // Otherwise calculate transmitted radiance
            return e + f.cwiseProduct(radiance(Ray(x,tdir), depth, Xi) * TP);
        }
    } else {
        // For depth <= 2, we compute the usual reflection + transmission
        return e + f.cwiseProduct(radiance(reflRay, depth, Xi) * Re) + 
                        f.cwiseProduct(radiance(Ray(x,tdir), depth, Xi) * Tr);
  }

}


#endif