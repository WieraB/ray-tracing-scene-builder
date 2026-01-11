#include <math.h> 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h> 
#include "../Eigen/Dense"
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

  // double intersect(const Ray &r) const { // returns distance, 0 if nohit
  //   Eigen::Vector3d op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
  //   double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
  //   if (det<0) return 0; else det=sqrt(det);
  //   return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  // }

  double intersect(const Ray &r, Eigen::Vector3d &n_out) const { // returns distance, 0 if nohit
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

struct Quadric_element {
  Eigen::Vector3d v0, v1, v2;
  double  A, B, C;  // Quadratic surface coefficients
  double D, E, F;
  double  G, H, I, J;
  Eigen::Vector3d e, c; // emission, color
  Refl_t refl;

  Quadric_element(double A_, double B_, double C_, 
    double D_, double E_, double F_, 
    double G_, double H_,  double I_, double J_,
    Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
    A(A_), B(B_), C(C_), 
    D(D_), E(E_), F(F_), 
    G(G_), H(H_), I(I_), J(J_),
    e(e_), c(c_), 
    refl(refl_) {}

  double intersect(const Ray &r, Eigen::Vector3d &n_out) const {
    const double eps = 1e-5;
    Eigen::Vector3d point;

    double  x0, y0, z0;
    double  dx, dy, dz;        
    double  A_q, B_q, C_q;        // Ray-quadric surface intersection coefficents
    double  dist;            // Distance to intersection
    double  root;            // Root of distance to intersection
    double t;

    dx = r.d[0];
    dy = r.d[1];
    dz = r.d[2];

    x0 = r.o[0];
    y0 = r.o[1];
    z0 = r.o[2];

    // F(x, y, z) = Ax^2 + 2Bxy + 2Cxz + 2Dx + Ey^2 + 2Fyz + 2Gy + 2Hz^2 + 2Iz + J = 0

    A_q = A * dx * dx + 2 * B * dx * dy + 2 * C * dx * dz +
                            E * dy * dy + 2 * F * dy * dz + 
                            H * dz * dz;
  
 
    B_q = 2 * (A * x0 * dx + B * (x0 * dy + dx * y0) + C * (x0 * dz + dx * z0) +
                        D * dx + E * y0 * dy + F * (y0 * dz + dy * z0) + G * dy +
                        H * z0 * dz + I * dz);
 
    C_q = A * x0 * x0 + 2 * B * x0 * y0 + 2 * C * x0 * z0 + 2 * D * x0 + 
                            E * y0 * y0 + 2 * F * y0 * z0 + 2 * G * y0 +
                            H * z0 * z0 + 2 * I * z0 + J;
 
   if ( A_q == 0.0 ) {  // It's a plane then.
      if ( B_q == 0.0 ) {
         return 0.0;
      }
 
      t = ( -C_q ) / ( B_q );
 
   } else {
      dist = B_q * B_q - 4 * A_q * C_q;
      if ( dist < 0.0 ) {
         return 0.0;
      }
 
      root = sqrt( dist );
      t = ( -B_q - root ) / ( A_q + A_q );
      if ( t < 0.0 ) {
         t = ( -B_q + root ) / ( A_q + A_q );
      }
   }
 
   if ( t < 0.001 )
      return 0.0;

   point = r.o + r.d * t;
   double x = point[0];
   double y = point[1];
   double z = point[2];

   n_out[0] = 2 * (A*x + B*y + C*z + D);
	 n_out[1] = 2 * (B*x + E*y + F*z + G);
	 n_out[2] = 2 * (C*x + F*y + H*z + I);
 
   return t;
  }
};

struct Quadratic_tet {
    std::vector<Eigen::Vector3d> nodes; // 10 points defining quadratic tet: 0-3 corners, 4-9 mid-edge nodes
    Eigen::Vector3d e, c; // emission, color
    Refl_t refl;

    Quadratic_tet(std::vector<Eigen::Vector3d> nodes_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_):
          nodes(nodes_), e(e_), c(c_), refl(refl_) {}

    // Quadratic shape functions for a tet
    static Eigen::VectorXd calc_shape_functions(double u, double v, double w) {
        double r = 1.0 - u - v - w;
        Eigen::VectorXd N(10);
        // For corner nodes
        N << r*(2*r-1), u*(2*u-1), v*(2*v-1), w*(2*w-1),
        // For mid-edge nodes
             4*u*r, 4*u*v, 4*v*r, 4*w*r, 4*u*w, 4*v*w;
        return N;
    }

    // Jacobians, i.e partial derivatives dP/du, dP/dv, dP/dw
    Eigen::Matrix3d calc_jacobians(double u, double v, double w) const {
        Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
        double r = 1.0 - u - v - w;

        // Derivatives of shape functions w.r.t u, v, w
        // This is a simplified representation of the gradient vectors
        double dNdu[10] = {1-4*r, 4*u-1, 0, 0, 4*(r-u), 4*v, -4*v, -4*w, 4*w, 0};
        double dNdv[10] = {1-4*r, 0, 4*v-1, 0, -4*u, 4*u, 4*(r-v), -4*w, 0, 4*w};
        double dNdw[10] = {1-4*r, 0, 0, 4*w-1, -4*u, 0, -4*v, 4*(r-w), 4*u, 4*v};

        for (int i = 0; i < 10; ++i) {
            Eigen::Vector3d pos = nodes[i];
            J.col(0) += pos * dNdu[i];
            J.col(1) += pos * dNdv[i];
            J.col(2) += pos * dNdw[i];
        }
        return J;
    }

    double intersect(const Ray &r, Eigen::Vector3d &n_out) const {

        // Initial guess is a tet's centroid
        Eigen::Vector3d uvw(0.25, 0.25, 0.25); 
        Eigen::Vector3d ro = r.o;
        Eigen::Vector3d rd = r.d;

        for (int iter = 0; iter < 10; ++iter) {
            Eigen::VectorXd N = calc_shape_functions(uvw.x(), uvw.y(), uvw.z());
            Eigen::Vector3d P = Eigen::Vector3d::Zero();
            for(int i=0; i<10; ++i) P += N[i] * nodes[i];

            // The objective is to minimise the distance between P(u,v,w) and Ray
            // We solve the system P(uvw) - (O + tD) = 0
            // Since we want t, we project the error onto the plane perpendicular to the ray
            Eigen::Matrix3d J = calc_jacobians(uvw.x(), uvw.y(), uvw.z());
            
            // To handle the ray, we solve for the intersection with the tangent plane
            Eigen::Vector3d residual = P - ro;
            double t = residual.dot(rd);
            Eigen::Vector3d f = residual - t * rd; // Perpendicular error

            if (f.norm() < 1e-7) {

              if (t > 0.0) {

                n_out[0] = uvw.x();
                n_out[1] = uvw.y();
                n_out[2] = uvw.z();
                return t;
              } else {
                return 0.0;
              }
            }

            // Newton step: J_reduced * delta_uvw = -f
            // We project J onto the 2D plane perpendicular to the ray
            Eigen::Matrix3d Proj = Eigen::Matrix3d::Identity() - rd * rd.transpose();
            Eigen::Matrix<double, 3, 3> reducedJ = Proj * J;
            
            // Solve using pseudo-inverse or SVD since it's an underdetermined system (3 vars to 2D plane)
            uvw -= reducedJ.completeOrthogonalDecomposition().solve(f);
        }

        return 0.0;

    }
    
};




// Sphere spheres[] = {//Scene: radius, position, emission, color, material
//   Sphere(600, Eigen::Vector3d(50, 681.6-.27, 81.6), Eigen::Vector3d(12,12,12), Eigen::Vector3d(0, 0, 0),     DIFF), //Light
//   Sphere(6, Eigen::Vector3d(-0.006, -5.0, -1.3),    Eigen::Vector3d(0, 0, 0),  Eigen::Vector3d(1,1,1)*.999,  REFR),//Glass
//   Sphere(100, Eigen::Vector3d(-200, -5.0, -200),    Eigen::Vector3d(0, 0, 0),  Eigen::Vector3d(.75,.75,.75), DIFF),//Back
// };

// Triangle triangles[] = {
//  Triangle(Eigen::Vector3d(20,20,40), Eigen::Vector3d(80,20,40), Eigen::Vector3d(50,60,60),
//           Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1,1,0)*.9, DIFF)
// };

// Quadratic_element quadratic_elements[] = {//Scene: radius, position, emission, color, material
//   Quadratic_element(1.0, 10.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 50.0, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };

// Quadric_element quadratic_elements[] = {//Scene: radius, position, emission, color, material
//   Quadric_element(1.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 1.0, 10, 50.0, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };

Quadric_element quadratic_elements[] = {//Scene: radius, position, emission, color, material
  
};



Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
Quadratic_tet({
  Eigen::Vector3d(0.0000000000000001,     -0.0000000000000000,      2.0000000000000000),
  Eigen::Vector3d(1.0271786786018278,     -1.4112658311988424,     -0.9763363743692252),
  Eigen::Vector3d(0.5468672504268858,     -1.9053669300092937,     -0.2655429766299571),
  Eigen::Vector3d(1.3388909735328707,     -1.4776893910037117,     -0.1542887705158133),
  Eigen::Vector3d(0.5135893393009140,     -0.7056329155994212,      0.5118318128153874),
  Eigen::Vector3d(0.2734336252134429,     -0.9526834650046468,      0.8672285116850215),
  Eigen::Vector3d(0.6694454867664354,     -0.7388446955018558,      0.9228556147420933),
  Eigen::Vector3d(0.8122941570524250,     -1.7115646775809932,     -0.6408779578637114),
  Eigen::Vector3d(1.2128643517229434,     -1.4808992624446611,     -0.5795665956648496),
  Eigen::Vector3d(0.9680931166509713,     -1.7367621660107437,     -0.2155293395495808)
}, 
  Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
};



std::vector<Sphere> spheres;

std::vector<Triangle> triangles;

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

inline bool intersect(const Ray &r, double &t, int &id, int &isTri, Eigen::Vector3d &n,
std::vector<Sphere> &spheres, std::vector<Triangle> triangles) {
  double inf=t=1e20;
  isTri=-1;
  double d;

  for (int i=0; i<(int)spheres.size(); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=spheres[i].intersect(r, n_tmp)) && d<t) { t=d; id=i; isTri=1; n=n_tmp; }
  }

  for (int i=0; i<(int)triangles.size(); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=triangles[i].intersect(r, n_tmp)) && d<t) {
      t=d; id=i; isTri=2; n=n_tmp;
    }
  }

  for (int i=0; i<(int)sizeof(quadratic_elements)/sizeof(Quadric_element); i++) {
    Eigen::Vector3d n_tmp;
    if ((d=quadratic_elements[i].intersect(r, n_tmp)) && d<t) {
      t=d; id=i; isTri=3; n=n_tmp;
      // std::cout << t; 
    }
  }

  for (int i=0; i<(int)sizeof(quadratic_tets)/sizeof(Quadratic_tet); i++) {
    // std::cout << "Hello";
    Eigen::Vector3d n_tmp;
    if ((d=quadratic_tets[i].intersect(r, n_tmp)) && d<t) {
      // std::cout << "Hello";
       t=d; id=i; isTri=4; n=n_tmp;
      
      Eigen::Vector3d n1 = n_tmp;

      if (n_tmp.allFinite()) {
          // n_tmp contains valid, usable values
          // std::cout << "yes";
          // std::cout << d;
          // std::cout << t;
          // std::cout << " "; 
      } else {
          // n_tmp has at least one NaN or Inf
          // std::cout << "no";
      }
   

      // std::cout << t; 
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

Refl_t stringToRefl(const std::string &s) {
    if (s == "DIFF" || s == "diff" || s == "diffuse") return DIFF;
    if (s == "SPEC" || s == "spec" || s == "specular") return SPEC;
    if (s == "REFR" || s == "refr" || s == "refract") return REFR;

    std::cerr << "Warning: unknown reflection type '" << s 
              << "', defaulting to DIFF\n";
    return DIFF;
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



Eigen::Vector3d radiance(const Ray &r, int depth, unsigned short *Xi){

  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  int obj_type=-1;
  Eigen::Vector3d n;

  if (!intersect(r, t, id, obj_type, n, spheres, triangles)) return Eigen::Vector3d(0, 0, 0); // if miss, return black
  Eigen::Vector3d x = r.o + r.d * t;
  Eigen::Vector3d nl = n.dot(r.d) < 0 ? n : n * -1; // Adjust normal direction
  Eigen::Vector3d f, e;
  Refl_t refl;

  if (obj_type==1) {
    const Sphere &obj=spheres[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  } else if (obj_type==2) {
    const Triangle &obj=triangles[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  } else if (obj_type==3) {
    const Quadric_element &obj=quadratic_elements[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  } else if (obj_type==4) {
    const Quadratic_tet &obj=quadratic_tets[id];
    f=obj.c; 
    e=obj.e; 
    refl=obj.refl;
  }

  // switch (obj_type)
  //   {
  //       case 1:
  //         f=spheres[id].c; 
  //         e=spheres[id].e; 
  //         refl=spheres[id].refl;
  //         std::cout << "type 1";
       
  //       case 2:        
  //         f=triangles[id].c; 
  //         e=triangles[id].e; 
  //         refl=triangles[id].refl;
  //       case 3:
  //         f=triangles[id].c; 
  //         e=triangles[id].e; 
  //         refl=triangles[id].refl;
  //   }
  
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
  int w=512,h=384,samps = argc >= 4 ? atoi (argv[3]) / 4 : 1;

  if(argc >= 2) loadOBJ(argv[1], triangles);
  if(argc >= 3) loadSpheres(argv[2], spheres);

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



  // std::vector<Eigen::Vector3d> nodes;

  // // p1      p2      p3      p4 etc
  // //1 25 43 38 227 228 229 163 153 162

  // nodes.push_back(Eigen::Vector3d(0.0000000000000001,     -0.0000000000000000,      2.0000000000000000));
  // nodes.push_back(Eigen::Vector3d(1.0271786786018278,     -1.4112658311988424,     -0.9763363743692252));
  // nodes.push_back(Eigen::Vector3d(0.5468672504268858,     -1.9053669300092937,     -0.2655429766299571));
  // nodes.push_back(Eigen::Vector3d(1.3388909735328707,     -1.4776893910037117,     -0.1542887705158133));
  // nodes.push_back(Eigen::Vector3d(0.5135893393009140,     -0.7056329155994212,      0.5118318128153874));
  // nodes.push_back(Eigen::Vector3d(0.2734336252134429,     -0.9526834650046468,      0.8672285116850215));
  // nodes.push_back(Eigen::Vector3d(0.6694454867664354,     -0.7388446955018558,      0.9228556147420933));
  // nodes.push_back(Eigen::Vector3d(0.8122941570524250,     -1.7115646775809932,     -0.6408779578637114));
  // nodes.push_back(Eigen::Vector3d(1.2128643517229434,     -1.4808992624446611,     -0.5795665956648496));
  // nodes.push_back(Eigen::Vector3d(0.9680931166509713,     -1.7367621660107437,     -0.2155293395495808));


  // Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
  // Quadratic_tet(nodes, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
  // };




}
