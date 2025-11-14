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
#include <memory>
#include <ranges>
#include <array>


#include "ray.h"
#include "aabb.h"

#define EIGEN_NO_DEBUG


// General base class for all objects in the scene
// Emission, colour, and reflection type should be specified for all objects.
// For the specific ray, intersect function calculates the distance to the object
// and the normal vector at the intersection point. 
class intersectable {

  public:
    Eigen::Vector3d e, c;      // emission, colour
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
    aabb bbox;  // bounding box
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
    intersectable(e_, c_, refl_), rad(rad_), p(p_) {
        Eigen::Vector3d disp = Eigen::Vector3d(rad_, rad_, rad_);
        bbox = aabb(p_ - disp, p_ + disp);
    }

    double intersect(const Ray &r, Eigen::Vector3d &n_out) const override { // returns distance, 0 if no hit
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
      intersectable(e_, c_, refl_), v0(v0_), v1(v1_), v2(v2_) {

        double x_min = std::min({v0_.x(), v1_.x(), v2_.x()});
        double x_max = std::max({v0_.x(), v1_.x(), v2_.x()});

        double y_min = std::min({v0_.y(), v1_.y(), v2_.y()});
        double y_max = std::max({v0_.y(), v1_.y(), v2_.y()});

        double z_min = std::min({v0_.z(), v1_.z(), v2_.z()});
        double z_max = std::max({v0_.z(), v1_.z(), v2_.z()});

        Eigen::Vector3d p_min = Eigen::Vector3d(x_min, y_min, z_min);
        Eigen::Vector3d p_max = Eigen::Vector3d(x_max, y_max, z_max);

        bbox = aabb(p_min, p_max);

      }
  
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

// List of objects in the scene
// Intersect function check if the specific ray intersects any objects in the scene,
// then it returns the information about the intersection closest to the camera.
class intersectable_list {
  private:
    aabb bbox;

  public:
  std::vector<std::shared_ptr<intersectable>> objects;
    intersectable_list() {};
    virtual ~intersectable_list() = default;

    void clear() { objects.clear(); }

    void add(std::shared_ptr<intersectable> object) {
      objects.push_back(object);
      bbox = aabb(bbox, object -> bbox);
    }

    int size() { return objects.size(); }

    const std::shared_ptr<intersectable>& operator[](size_t i) const {
        return objects[i];
    }

    inline bool intersect(const Ray &r, double &t, int &id, Eigen::Vector3d &n) {
      double inf=t=1e20;
      double d;
    
      for (int i=0; i<(int)size(); i++) {
        Eigen::Vector3d n_tmp;
        if ((d=objects[i] -> intersect(r, n_tmp)) && d<t) {
           t=d; id=i; n=n_tmp;
        }
      }
    
      return t<inf;
    }

    // void add_block(std::vector<std::shared_ptr<intersectable>> objects_) {
    //   objects.insert( objects.end(), objects_.begin(), objects_.end() );
    // }

};


#endif