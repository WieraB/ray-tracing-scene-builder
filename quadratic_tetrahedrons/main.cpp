#include <math.h> 
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h> 
#include "../Eigen/Dense"
// #include "../Eigen/Eigenvalues" 
#include "dlib/optimization.h"
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
    std::vector<Eigen::Vector3d> nodes; 
    Eigen::Vector3d e, c;
    Refl_t refl;

    // Nodes for the 4 faces - each face has 6 nodes: 3 corners, 3 mid-edges
    // Numbering for 1 tetrahedron example
    // const int face_indices[4][6] = {
    //     {0, 2, 1, 6, 5, 4}, // Face 0 (bottom)
    //     {0, 1, 3, 4, 8, 7}, // Face 1
    //     {1, 2, 3, 5, 9, 8}, // Face 2
    //     {2, 0, 3, 6, 7, 9}  // Face 3
    // };

      // Correct numbering for NGSolve/Netgen
      const int face_indices[4][6] = {
      {0, 1, 2, 4, 7, 5}, // Face 0 (bottom) (v0, v1, v2)
      {0, 1, 3, 4, 8, 6}, // Face 1 (v0, v1, v3)
      {1, 2, 3, 7, 9, 8}, // Face 2 (v1, v2, v3)
      {0, 2, 3, 5, 9, 6}  // Face 3 (v0, v2, v3)
      };

    Quadratic_tet(std::vector<Eigen::Vector3d> nodes_, Eigen::Vector3d e_, Eigen::Vector3d c_, Refl_t refl_) :
        nodes(nodes_), e(e_), c(c_), refl(refl_) {}

        // Quadratic triangle shape functions (g, h)
        // r = 1 - g - h
        static Eigen::VectorXd get_face_N(double g, double h) {
        double r = 1.0 - g - h;
        Eigen::VectorXd N(6);
        N << r*(2*r-1), g*(2*g-1), h*(2*h-1), 4*g*r, 4*g*h, 4*h*r;
        return N;
    }

    static Eigen::Matrix<double, 3, 2> get_face_Jacobian(double g, double h, const std::vector<Eigen::Vector3d>& f_nodes) {
        double r = 1.0 - g - h;
        // Derivatives of N wrt g and h
        double dNdu[6] = { -(4*r-1), 4*g-1, 0, 4*(r-g), 4*h, -4*h };
        double dNdv[6] = { -(4*r-1), 0, 4*h-1, -4*g, 4*g, 4*(r-h) };

        Eigen::Matrix<double, 3, 2> J = Eigen::Matrix<double, 3, 2>::Zero();
        for (int i = 0; i < 6; ++i) {
            J.col(0) += f_nodes[i] * dNdu[i];
            J.col(1) += f_nodes[i] * dNdv[i];
        }
        return J;
    }

    static Eigen::Matrix<double, 3, 3> get_face_Hessian(const std::vector<Eigen::Vector3d>& f_nodes) {
        // Constant second derivatives of the 6 shape functions
        // d2N/dg2
        double ddN_dg2[6] = { 4.0, 4.0, 0.0, -8.0, 0.0, 0.0 };
        // d2N/dh2
        double ddN_dh2[6] = { 4.0, 0.0, 4.0, 0.0, 0.0, -8.0 };
        // d2N/dgdh
        double ddN_dgh[6] = { 4.0, 0.0, 0.0, -4.0, 4.0, -4.0 };
    
        Eigen::Matrix<double, 3, 3> H = Eigen::Matrix<double, 3, 3>::Zero();
        Eigen::Vector3d d2N_dg2 = Eigen::Vector3d::Zero();
        Eigen::Vector3d d2N_dh2 = Eigen::Vector3d::Zero();
        Eigen::Vector3d d2N_dgh = Eigen::Vector3d::Zero();
    
        for (int i = 0; i < 6; ++i) {
            d2N_dg2 += f_nodes[i] * ddN_dg2[i];
            d2N_dh2 += f_nodes[i] * ddN_dh2[i];
            d2N_dgh += f_nodes[i] * ddN_dgh[i];
        }

        H.col(0) = d2N_dg2; // d2N/dg2
        H.col(1) = d2N_dh2; // d2N/dh2
        H.col(2) = d2N_dgh; // d2N/dgdh
        
        return H;
    }


    static Eigen::Vector2d get_eq_Jacobian(double g, double h, const std::vector<Eigen::Vector3d>& f_nodes,
      const Ray &ray) {
        
        // Calculate the position on the surface P
        Eigen::VectorXd N = get_face_N(g, h);
        Eigen::Vector3d P = Eigen::Vector3d::Zero();
        for(int i=0; i<6; ++i) P += N[i] * f_nodes[i];

        // Calculate Jacobians of the tet face shape functions
        Eigen::Matrix<double, 3, 2> J = get_face_Jacobian(g, h, f_nodes);

        // Calculate V and t
        Eigen::Vector3d V = P - ray.o;
        double t = V.dot(ray.d);

        // Derivatives of F wrt g and h
        double dFdg = 2 * J.col(0).dot(V - t * ray.d);
        double dFdh = 2 * J.col(1).dot(V - t * ray.d);

        Eigen::Vector2d J_F(dFdg, dFdh);

        return J_F;
    }


    static Eigen::Matrix<double, 2, 2> get_eq_Hessian(double g, double h, const std::vector<Eigen::Vector3d>& f_nodes,
      const Ray &ray) {
        
        // Calculate the position on the surface P
        Eigen::VectorXd N = get_face_N(g, h);
        Eigen::Vector3d P = Eigen::Vector3d::Zero();
        for(int i=0; i<6; ++i) P += N[i] * f_nodes[i];

        // Calculate Jacobians of the tet face shape functions
        Eigen::Matrix<double, 3, 2> J = get_face_Jacobian(g, h, f_nodes);

        // Calculate Jacobians of the tet face shape functions
        Eigen::Matrix<double, 3, 3> H = get_face_Hessian(f_nodes);

        // Calculate V and t
        Eigen::Vector3d V = P - ray.o;
        double t = V.dot(ray.d);
        Eigen::Vector3d R = ray.o + t * ray.d;

        Eigen::Vector3d unity = Eigen::Vector3d::Ones();

        // Derivatives of F wrt g and h
        // double ddF_dg2 = 2 * H.col(0).dot(P + unity - R);
        // double ddF_dh2 = 2 * H.col(1).dot(P + unity - R);
        // double ddF_dgh = 2 * H.col(2).dot(P + unity - R);

        // double ddF_dg2 = 2 * H.col(0).dot(V + unity) * (1 - ray.d.dot(ray.d));
        // double ddF_dh2 = 2 * H.col(1).dot(V + unity) * (1 - ray.d.dot(ray.d));
        // double ddF_dgh = 2 * H.col(2).dot(V + unity) * (1 - ray.d.dot(ray.d));

        // double ddF_dg2 = 2 * (H.col(0).dot(V) + J.col(0).dot(J.col(0)))
        //                - 2 * ray.d.dot(ray.d) * (H.col(0).dot(V) + J.col(0).dot(J.col(0)));

        // double ddF_dh2 = 2 * (H.col(1).dot(V) + J.col(1).dot(J.col(1)))
        //                - 2 * ray.d.dot(ray.d) * (H.col(1).dot(V) + J.col(1).dot(J.col(1)));

        // double ddF_dgh = 2 * (H.col(2).dot(V) + J.col(0).dot(J.col(1)))
        //                - 2 * ray.d.dot(ray.d) * (H.col(2).dot(V) + J.col(0).dot(J.col(1)));


        double ddF_dg2 = 2 * (J.col(0).dot(J.col(0)) - (J.col(0).dot(ray.d)) * (J.col(0).dot(ray.d)) 
                            + V.dot(H.col(0)) - V.dot(ray.d) * H.col(0).dot(ray.d));

        double ddF_dh2 = 2 * (J.col(1).dot(J.col(1)) - (J.col(1).dot(ray.d)) * (J.col(1).dot(ray.d)) 
                            + V.dot(H.col(1)) - V.dot(ray.d) * H.col(1).dot(ray.d));

        double ddF_dgh = 2 * (J.col(0).dot(J.col(1)) - (J.col(0).dot(ray.d)) * (J.col(1).dot(ray.d)) 
                            + V.dot(H.col(2)) - V.dot(ray.d) * H.col(2).dot(ray.d));

        Eigen::Matrix<double, 2, 2> H_F 
        {
          {ddF_dg2, ddF_dgh},
          {ddF_dgh, ddF_dh2},
        };

        return H_F;
    }


    double backtracking_line_search(Eigen::Vector2d dir, Eigen::Vector2d x_init, Eigen::Vector2d J_F,
      const std::vector<Eigen::Vector3d>& f_nodes, const Ray &r) const {

      double s = 1.0;          // Initial step size
      // double beta = 0.5;       // Contraction factor (0 < beta < 1)
      // double c = 1e-4;         // Armijo condition constant
      // int max_iter = 100;

      double beta = 0.5;       // Contraction factor (0 < beta < 1)
      double c = 1e-4;         // Armijo condition constant
      int max_iter = 1000;

      Eigen::VectorXd N = get_face_N(x_init.x(), x_init.y());
      Eigen::Vector3d P = Eigen::Vector3d::Zero();
      for(int i=0; i<6; ++i) P += N[i] * f_nodes[i];
      Eigen::Vector3d V = P - r.o;
      double dotProduct = V.dot(r.d);
      double fx = V.squaredNorm() - (dotProduct * dotProduct);



      double gradient_dir = J_F.dot(dir); // Directional derivative

      //If the direction isn't a descent direction, line search will fail
      if (gradient_dir > 0) return 0.0;


      Eigen::Vector2d x = x_init;
      for (int i = 0; i < max_iter; ++i) {
        Eigen::Vector2d x_next = x + s * dir;

        N = get_face_N(x_next.x(), x_next.y());
        P = Eigen::Vector3d::Zero();
        for(int i=0; i<6; ++i) P += N[i] * f_nodes[i];
        V = P - r.o;
        dotProduct = V.dot(r.d);
        double f_next = V.squaredNorm() - (dotProduct * dotProduct);

        // Armijo Condition: F(x + s*p) <= F(x) + c * s * (grad_F^T * p)
        if (f_next <= fx + c * s * gradient_dir) {
            break; // Success: sufficient decrease found
        }

        s *= beta; // Shrink the step size

        if (s < 1e-12) break; // Step size too small
    }

    x = x + s * dir;
    // std::cout << "Optimal scalar s: " << s << std::endl;
    // std::cout << "New x: \n" << (x).transpose() << std::endl;

    return s;

    }

    double intersect(const Ray &r, Eigen::Vector3d &n_out) const {
      
      // Define imprecision parameters for the initial guess, they are especially needed for the cases
      // where the ray is close to being tangential to the linear triangle.
      // If the determinant is close to zero, the ray misses the triangle.
      // Imprecision levels should be relatively large for the initial guess.
      const double eps_init_guess1 = 4e-1; // Imprecision for linear tringle's determinant
      const double eps_init_guess2 = 4e-1; // Imprecision for linear tringle's isoparametric coordinates (u, v)

      const double eps_t = 1e-8; // Imprecision for t

      const int iter_max = 1500; // Maximum number of iterations for Newton-Raphson method.

      // Define imprecision parameters for the solution.
      // Imprecision levels should be relatively small for solution.
      const double eps_sol1 = 1; // Imprecision for distance
      const double eps_sol2 = 1e-8; // Imprecision for the quadratic triangle's isoparametric coordinates (g, h)
      const double eps_sol3 = 1e-10; // Imprecision for the determinant
      const double eps_sol4 = 0.1; // Convergence criteria, max allowed relative difference


      double min_t = 1e20;
      bool intersect = false;

        // Iterate through each of the 4 faces
        for (int f = 0; f < 4; ++f) {
            std::vector<Eigen::Vector3d> f_nodes(6);
            for(int i=0; i<6; ++i) f_nodes[i] = nodes[face_indices[f][i]];

            // 1. Intersect ray with the plane defined by 3 corner nodes to get the initial guess for t
            // Use the code already used for Traingle structure.
            // Impose the defined imprecision levels
            Eigen::Vector3d v0 = f_nodes[0], v1 = f_nodes[1], v2 = f_nodes[2];
            Eigen::Vector3d edge1 = v1 - v0, edge2 = v2 - v0;
            Eigen::Vector3d pvec = r.d.cross(edge2);
            double det = edge1.dot(pvec);
            if (fabs(det) < eps_init_guess1) continue;
            double invDet = 1.0 / det;
            Eigen::Vector3d tvec = r.o - v0;
            double u = tvec.dot(pvec) * invDet;
            if (u < 0 - eps_init_guess2 || u > 1 + eps_init_guess2) continue;
            Eigen::Vector3d qvec = tvec.cross(edge1);
            double v = r.d.dot(qvec) * invDet;
            if (v < 0 - eps_init_guess2 || u + v > 1 + eps_init_guess2) continue;
            double t_guess = edge2.dot(qvec) * invDet;
            if (t_guess < eps_t) continue;

            // 2. Start Newton-Raphson at t_guess, and gh from the right triangle's centroid 
            // or approximated (u, v) from the linear triangle. 
            // Eigen::Vector2d gh(0.33, 0.33); 
            Eigen::Vector2d gh(u, v); 
            double t = t_guess;
            Eigen::Matrix<double, 2, 2> M = get_eq_Hessian(gh.x(), gh.y(), f_nodes, r);

            double res = 1e20;
            double res_prev = 1e20;

            for (int iter = 0; iter < iter_max; ++iter) {
                Eigen::VectorXd N = get_face_N(gh.x(), gh.y());
                Eigen::Vector3d P = Eigen::Vector3d::Zero();
                for(int i=0; i<6; ++i) P += N[i] * f_nodes[i];

                Eigen::Vector3d V = P - r.o;
                t = V.dot(r.d);
                res_prev = res;
                res = V.squaredNorm() - (t * t);
                double diff = abs(res-res_prev) / abs(res);

                if ( diff < eps_sol4) {

                    if (res >= eps_sol1) break;
                    // std::cout << diff << "    ";
                  
                    if (gh.x() >= 0 - eps_sol2 && gh.y() >= 0 - eps_sol2 && (gh.x() + gh.y()) <= 1 + eps_sol2) {
                        if (t < min_t && t > eps_t) {
                            min_t = t;
                            Eigen::Matrix<double, 3, 2> J = get_face_Jacobian(gh.x(), gh.y(), f_nodes);
                            n_out = J.col(0).cross(J.col(1)).normalized();
                            intersect = true;
                        }
                    }
                    break;
                }

                // Solve H_F * [dg, dh]^T = -J_F
                M = get_eq_Hessian(gh.x(), gh.y(), f_nodes, r);
                Eigen::Vector2d J_F = get_eq_Jacobian(gh.x(), gh.y(), f_nodes, r);
                Eigen::Vector2d dir = M.inverse() * (-J_F);

                double s = backtracking_line_search(dir, gh, J_F, f_nodes, r);

                if (std::abs(M.determinant()) < eps_sol3) break;
                Eigen::Vector2d delta = s * M.inverse() * (-J_F);
                gh.x() += delta.x();
                gh.y() += delta.y();
            }
        }

        return intersect ? min_t : 0.0;
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

// Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
  
// };



// Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
// Quadratic_tet({
//   Eigen::Vector3d(0.0000000000000001,     -0.0000000000000000,      2.0000000000000000),
//   Eigen::Vector3d(1.0271786786018278,     -1.4112658311988424,     -0.9763363743692252),
//   Eigen::Vector3d(0.5468672504268858,     -1.9053669300092937,     -0.2655429766299571),
//   Eigen::Vector3d(1.3388909735328707,     -1.4776893910037117,     -0.1542887705158133),
//   Eigen::Vector3d(0.5135893393009140,     -0.7056329155994212,      0.5118318128153874),
//   Eigen::Vector3d(0.2734336252134429,     -0.9526834650046468,      0.8672285116850215),
//   Eigen::Vector3d(0.6694454867664354,     -0.7388446955018558,      0.9228556147420933),
//   Eigen::Vector3d(0.8122941570524250,     -1.7115646775809932,     -0.6408779578637114),
//   Eigen::Vector3d(1.2128643517229434,     -1.4808992624446611,     -0.5795665956648496),
//   Eigen::Vector3d(0.9680931166509713,     -1.7367621660107437,     -0.2155293395495808)
// }, 
//   Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };


// Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
// Quadratic_tet({
//   Eigen::Vector3d(0.943,	0.000,	-0.333),
//   Eigen::Vector3d(-0.471,	0.817,	-0.333),
//   Eigen::Vector3d(-0.471,	-0.817,	-0.333),
//   Eigen::Vector3d(0.000,	0.000,	1.000),
//   Eigen::Vector3d(0.236,	0.408,	-0.333),
//   Eigen::Vector3d(-0.471,	0.000,	-0.333),
//   Eigen::Vector3d(0.236,	-0.408,	-0.333),
//   Eigen::Vector3d(0.471,	0.000,	0.333),
//   Eigen::Vector3d(-0.236,	0.408,	0.333),
//   Eigen::Vector3d(-0.236,	-0.408,	0.333)
// }, 
//   Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };



// Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
// Quadratic_tet({
//   Eigen::Vector3d(6.000,	0.000,	-2.121),
//   Eigen::Vector3d(-3.000,	5.196,	-2.121),
//   Eigen::Vector3d(-3.000,	-5.196,	-2.121),
//   Eigen::Vector3d(0.000,	0.000,	6.000),
//   Eigen::Vector3d(1.500,	2.598,	-2.121),
//   Eigen::Vector3d(-3.000,	0.000,	-2.121),
//   Eigen::Vector3d(1.500,	-2.598,	-2.121),
//   Eigen::Vector3d(3.000,	0.000,	1.940),
//   Eigen::Vector3d(-1.5,	2.598,	1.940),
//   Eigen::Vector3d(-1.5,	-2.598,	1.940)
// }, 
//   Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };


// Quadratic_tet quadratic_tets[] = {//Scene: radius, position, emission, color, material
// Quadratic_tet({
//   Eigen::Vector3d(6.000,	0.000,	-2.121),
//   Eigen::Vector3d(-3.000,	5.196,	-2.121),
//   Eigen::Vector3d(-3.000,	-5.196,	-2.121),
//   Eigen::Vector3d(0.000,	0.000,	6.000),
//   Eigen::Vector3d(1.500,	2.598,	-2.121) + Eigen::Vector3d(0.00,	0.00,	2.0),
//   Eigen::Vector3d(-3.000,	0.000,	-2.121),
//   Eigen::Vector3d(1.500,	-2.598,	-2.121) + Eigen::Vector3d(0.00,	0.00,	2.0),
//   Eigen::Vector3d(3.000,	0.000,	1.940),
//   Eigen::Vector3d(-1.5,	2.598,	1.940),
//   Eigen::Vector3d(-1.5,	-2.598,	1.940)
// }, 
//   Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF)
// };







std::vector<Sphere> spheres;

std::vector<Triangle> triangles;

std::vector<Quadratic_tet> quadratic_tets;

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

  for (int i=0; i<(int)quadratic_tets.size(); i++) {
    // std::cout << "HELLO";
    Eigen::Vector3d n_tmp;
    if ((d=quadratic_tets[i].intersect(r, n_tmp)) && d<t) {
       t=d; id=i; isTri=4; n=n_tmp;
    }
  }

  return t<inf;
}


struct Mesh {
    std::vector<Eigen::Vector3d> points;
    std::vector<std::vector<int>> elements;
    std::vector<std::vector<Eigen::Vector3d>> elem_coords;


    void loadVolFile(const std::string& filename) {
      // Mesh mesh;
      std::ifstream file(filename);
      
      if (!file.is_open()) {
          std::cerr << "Cannot open .vol file: " << filename << std::endl;
      }
  
      std::string line;
      while (std::getline(file, line)) {
          // Load points
          if (line.find("points") != std::string::npos) {
            // std::cout << "yes1  ";
              int numPoints;
              file >> numPoints;
              points.reserve(numPoints);
              for (int i = 0; i < numPoints; ++i) {
                  double x, y, z;
                  file >> x >> y >> z;
                  points.emplace_back(x, y, z);
              }
  
              break;
          }
      }
      file.close();
  
  
      std::ifstream file1(filename);
      std::string line1;
      while (std::getline(file1, line1)) {
          
          // Load elements
          if (line1.find("volumeelements") != std::string::npos) {
            // std::cout << "yes2!   ";
              int numElements;
              file1 >> numElements;
              elements.reserve(numElements);
              for (int i = 0; i < numElements; ++i) {
                  int mat, np;
                  file1 >> mat >> np; // material number of points
                  
                  std::vector<int> elementNodes(np);
                  for (int j = 0; j < np; ++j) {
                      file1 >> elementNodes[j];
                      // Convert 1-based indexing used in .vol to 0-based
                      elementNodes[j]--; 
                  }
                  elements.push_back(elementNodes);

                  // if (i == 1) break;
              }
          }
      }
      file1.close();
    };

    void getElementCoords() {

      std::vector<int> element_node_numbers;

      elem_coords.resize(elements.size()); // Number of rows
      for(auto& row : elem_coords) {
          row.resize(10, Eigen::Vector3d::Zero()); // Number of columns
      }

      for (int i = 0; i < elements.size(); ++i) {

        element_node_numbers = elements[i];


        for (int j = 0; j < 10; ++j) {

          // elem_coords[i][j] = Eigen::Vector3d(0, 0, 1);
          elem_coords[i][j] = points[element_node_numbers[j]];
        }

      }

    };
};


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

bool loadVOL(const Mesh &mesh , std::vector<Quadratic_tet> &quadratic_tets) {
  
  std::vector<Eigen::Vector3d> nodes; 

  for (int i = 0; i < mesh.elem_coords.size(); ++i) {

    nodes = mesh.elem_coords[i];
    quadratic_tets.emplace_back(nodes, Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(1.0, 0.078, 0.576), DIFF);

  }
  
  std::cerr << "Loaded " << quadratic_tets.size() << " quadratic tetrahedrons" << "\n";
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



   Mesh mesh;
   mesh.loadVolFile("sphere.vol");
   mesh.getElementCoords();



   std::cout << "Number of points: " << mesh.points.size() << std::endl;
   std::cout << "Number of elements: " << mesh.elements.size() << std::endl;
   std::cout << "Number of elements: " << mesh.elem_coords.size() << std::endl;

   loadVOL(mesh , quadratic_tets);

   std::cout << "Number of loaded elements: " << quadratic_tets.size() << std::endl;





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
