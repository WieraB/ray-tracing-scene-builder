#ifndef LOAD_H
#define LOAD_H

#include <math.h> 
#include <stdlib.h>
#include <stdio.h> 
#include "../Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "hittable.h"

#define EIGEN_NO_DEBUG

bool loadOBJ(const std::string &filename, std::vector<Triangle> &tris) {
  std::ifstream in(filename);
  if(!in) {
    std::cerr << "Cannot open OBJ file: " << filename << "\n";
    return false;
  }
  std::vector<Eigen::Vector3d> verts;
  std::string line;
  while(std::getline(in, line)) {
    if(line.size()<2) continue;
    std::istringstream iss(line);
    if(line[0]=='v' && line[1]==' ') {
      double x,y,z; iss.ignore(2); iss >> x >> y >> z;
      verts.emplace_back(x,y,z);
    }
    else if(line[0]=='f') {
      std::string v1,v2,v3;
      iss.ignore(2); iss >> v1 >> v2 >> v3;
      auto idx = [&](const std::string &s) {
        return atoi(s.substr(0, s.find('/')).c_str()) - 1;
      };
      int i0=idx(v1), i1=idx(v2), i2=idx(v3);
      if(i0<0||i1<0||i2<0||i0>=verts.size()||i1>=verts.size()||i2>=verts.size()) continue;
      tris.emplace_back(verts[i0], verts[i1], verts[i2],
                        Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1,1,0)*.8, DIFF);
    }
  }
  std::cerr << "Loaded " << tris.size() << " triangles from " << filename << "\n";
  return true;
}

bool loadSpheres(const std::string &filename, std::vector<Sphere> &sph) {
  std::ifstream in(filename);
  if(!in) {
    std::cerr << "Cannot open OBJ file: " << filename << "\n";
    return false;
  }
  std::vector<Eigen::Vector3d> verts;
  std::string line;
  while(std::getline(in, line)) {
    if(line.size()<2) continue;
    std::istringstream iss(line);
    if(line[0]=='s' && line[1]==' ') {
      double rad, px, py, pz, ex, ey, ez, cx, cy, cz;
      std::string refl_s;

      iss.ignore(2); 
      iss >> rad;
      iss >> px >> py >> pz;
      iss >> ex >> ey >> ez;
      iss >> cx >> cy >> cz;
      iss >> refl_s;

      sph.emplace_back(rad, Eigen::Vector3d(px, py, pz), Eigen::Vector3d(ex, ey, ez), Eigen::Vector3d(cx, cy, cz), stringToRefl(refl_s));
    }
  }
  
  std::cerr << "Loaded " << sph.size() << " spheres from " << filename << "\n";
  return true;
}

#endif