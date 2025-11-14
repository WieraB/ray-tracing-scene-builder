#ifndef RAY_H
#define RAY_H

#include <math.h> 
#include <stdlib.h>
#include <stdio.h> 
#include "../Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#define EIGEN_NO_DEBUG

struct Ray { Eigen::Vector3d o, d; Ray(Eigen::Vector3d o_, Eigen::Vector3d d_) : o(o_), d(d_) {} };

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

Refl_t stringToRefl(const std::string &s) {
    if (s == "DIFF" || s == "diff" || s == "diffuse") return DIFF;
    if (s == "SPEC" || s == "spec" || s == "specular") return SPEC;
    if (s == "REFR" || s == "refr" || s == "refract") return REFR;

    std::cerr << "Warning: unknown reflection type '" << s 
              << "', defaulting to DIFF\n";
    return DIFF;
}


#endif