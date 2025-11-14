#include <math.h> 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h> 
#include "../Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

#include "intersectable.h"
#include "ray.h"
#include "load.h"
#include "misc.h"
#include "camera.h"
#include "aabb.h"
#include "bvh.h"

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
#define EIGEN_NO_DEBUG

int main(int argc,char *argv[]){
  int w=512,h=384,samps = argc >= 4 ? atoi (argv[3]) / 4 : 1;

  intersectable_list scene;

  if(argc >= 2) loadOBJ(argv[1], scene);
  if(argc >= 3) loadSpheres(argv[2], scene);

  // auto world = intersectable_list(std::make_shared<bvh_node>(scene));

  // scene = intersectable_list(std::make_shared<bvh_node>(scene));

  std::cout << "Scene size : " << scene.size();

  double fov = 20;  // Vertical view angle (field of view)
  auto theta = degrees_to_radians(fov);
  auto scale = std::tan(theta/2);
  double aspect_ratio = double(w) / double(h);

  Eigen::Vector3d camPos(110, -10, 110);
  Eigen::Vector3d camDir = (Eigen::Vector3d(0, 0, 0) - camPos).normalized(); // look at origin

  Ray cam(camPos, camDir);

  Eigen::Vector3d right = cam.d.cross(Eigen::Vector3d(0, 1, 0)).normalized(); // orthonormal camera basis
  Eigen::Vector3d up = right.cross(cam.d).normalized(); // orthonormal camera basis

  Eigen::Vector3d cx = right * aspect_ratio * scale; // scale the orthogonal basis
  Eigen::Vector3d cy = up * scale; // scale the orthogonal basis
  
  Eigen::Vector3d r=Eigen::Vector3d(0, 0, 0);
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> c(w*h);

  unsigned short yyy;
  #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y=0; y<h; y++){                       // Loop over image rows
      fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
      yyy = y*y*y;
      for (unsigned short x=0, Xi[3]={0,0,yyy}; x<w; x++)   // Loop cols
        for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
          for (int sx=0; sx<2; sx++, r=Eigen::Vector3d(0, 0, 0)){        // 2x2 subpixel cols
            
            for (int s=0; s<samps; s++){

              double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
              double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
              Eigen::Vector3d d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                      cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
              r = r + radiance(Ray(cam.o+d*140,d.normalized()),0,Xi, scene)*(1./samps);
            }
            c[i] = c[i] + Eigen::Vector3d(clamp(r.x()),clamp(r.y()),clamp(r.z()))*.25;
          }
  }

  FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(f,"%d %d %d ", toInt(c[i].x()), toInt(c[i].y()), toInt(c[i].z()));
  fclose(f);
  fprintf(stderr,"\nDone. Wrote image.ppm\n");

}
