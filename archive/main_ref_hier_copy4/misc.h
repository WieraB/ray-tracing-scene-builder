#ifndef MISC_H
#define MISC_H

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

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

#endif