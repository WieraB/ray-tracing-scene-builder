#include <math.h> 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h> 
#include "./Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

// #define EIGEN_USE_MKL_ALL
// #define EIGEN_USE_BLAS
#define EIGEN_NO_DEBUG

// Add trinagles to the schene

struct Ray { Eigen::Vector3d o, d; Ray(Eigen::Vector3d o_, Eigen::Vector3d d_) : o(o_), d(d_) {} };

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere {
  double rad;       // radius

  Eigen::Vector3d p, e, c;      // position, emission, color

  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

  Sphere(double rad_, Eigen::Vector3d p_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Eigen::Vector3d op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
    if (det<0) return 0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};

struct Triangle {
  Eigen::Vector3d v0, v1, v2;
  Eigen::Vector3d e, c; // emission, color
  Refl_t refl;

  Triangle(Eigen::Vector3d v0_, Eigen::Vector3d v1_, Eigen::Vector3d v2_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
    v0(v0_), v1(v1_), v2(v2_), e(e_), c(c_), refl(refl_) {}

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

Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(600, Eigen::Vector3d(50, 681.6-.27, 81.6),Eigen::Vector3d(12,12,12), Eigen::Vector3d(0, 0, 0),     DIFF), //Light
  Sphere(6, Eigen::Vector3d(-0.006, -5.0, -1.3),    Eigen::Vector3d(0, 0, 0),Eigen::Vector3d(1,1,1)*.999, REFR),//Glass
  Sphere(100, Eigen::Vector3d(-200, -5.0, -200),     Eigen::Vector3d(0, 0, 0),  Eigen::Vector3d(.75,.75,.75), DIFF),//Back
};

// Triangle triangles[] = {
//  Triangle(Eigen::Vector3d(20,20,40), Eigen::Vector3d(80,20,40), Eigen::Vector3d(50,60,60),
//           Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1,1,0)*.9, DIFF)
// };

std::vector<Triangle> triangles;

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

inline bool intersect(const Ray &r, double &t, int &id, bool &isTri, Eigen::Vector3d &nTri) {
  double inf=t=1e20;
  isTri=false;
  double d;
  for (int i=0; i<sizeof(spheres)/sizeof(Sphere); i++) {
    if ((d=spheres[i].intersect(r)) && d<t) { t=d; id=i; isTri=false; }
  }
  for (int i=0; i<(int)triangles.size(); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=triangles[i].intersect(r, n_tmp)) && d<t) {
      t=d; id=i; isTri=true; nTri=n_tmp;
    }
  }
  return t<inf;
}


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


Eigen::Vector3d radiance(const Ray &r, int depth, unsigned short *Xi){

  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  bool isTri=false;
  Eigen::Vector3d nTri;
  if (!intersect(r, t, id, isTri, nTri)) return Eigen::Vector3d(0, 0, 0); // if miss, return black
  Eigen::Vector3d x = r.o + r.d * t, n, nl, f, e;
  Refl_t refl;

  if (!isTri) {
    const Sphere &obj=spheres[id];
    n = (x - obj.p).normalized();    // Normal at the intersection point
    nl = n.dot(r.d) < 0 ? n : n * -1; // Adjust normal direction
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  } else {
    const Triangle &tri=triangles[id];
    n=nTri; 
    nl=n.dot(r.d) < 0 ? n : n * -1;
    f=tri.c; e=tri.e; refl=tri.refl;
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

int main(int argc,char *argv[]){
  int w=512,h=384,samps = argc >= 3 ? atoi (argv[2]) / 4 : 1;

  if(argc >= 2) loadOBJ(argv[1], triangles);

  double fov = 20;  // Vertical view angle (field of view)
  auto theta = degrees_to_radians(fov);
  auto scale = std::tan(theta/2);
  double aspect_ratio = double(w) / double(h);

  Eigen::Vector3d camPos(110, -10, 110);
  Eigen::Vector3d camDir = (Eigen::Vector3d(0, 0, 0) - camPos).normalized(); // look at origin

  Ray cam(camPos, camDir);

  Eigen::Vector3d right = cam.d.cross(Eigen::Vector3d(0, 1, 0)).normalized();
  Eigen::Vector3d up = right.cross(cam.d).normalized();

  Eigen::Vector3d cx = right * aspect_ratio * scale;
  Eigen::Vector3d cy = up * scale;
  
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
              r = r + radiance(Ray(cam.o+d*140,d.normalized()),0,Xi)*(1./samps);
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
