#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2


// The `Vec` struct represents a 3D vector and is used for positions, directions, and colors. It provides several operations for vector math.
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
  double x, y, z;                  // position, also color (r,g,b)
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }

  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }

  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }

  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }

  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }

  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }

  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:

  Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};

//  The `Ray` struct represents a ray, which has an origin (`o`) and a direction (`d`). A ray is used to trace light paths in the scene.
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

// This enum defines the types of materials:
// `DIFF` (Diffuse) for matte surfaces.
// `SPEC` (Specular) for reflective surfaces.
// `REFR` (Refraction) for transparent surfaces (e.g., glass).
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

// The `Sphere` struct defines the properties of a sphere: 
// its radius (`rad`), position (`p`), emission color (`e`), surface color (`c`), and reflection type (`refl`).
// The `intersect` function computes the intersection of a ray with the sphere. 
// It solves the quadratic equation for ray-sphere intersection and returns the distance to the intersection point.
struct Sphere {
  double rad;       // radius

  Vec p, e, c;      // position, emission, color

  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
    if (det<0) return 0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};

// This is an array of `Sphere` objects that represents the scene. 
// Each sphere is defined by its radius, position, emission color, surface color, and reflection type. 
// These spheres make up the environment the light rays will interact with.
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};

// `clamp`: Clamps a value between 0 and 1.
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

// `toInt`: Converts a normalized value (between 0 and 1) to an integer in the range [0, 255] for image colors.
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

// The `intersect` function checks for intersections between a ray and the scene (represented by the `spheres` array). 
// It finds the closest sphere the ray intersects and returns the distance (`t`) and the ID of the intersected sphere.
inline bool intersect(const Ray &r, double &t, int &id){
  double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
  for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
  return t<inf;
}

// The `radiance` function recursively calculates the colour of a pixel based on how rays interact with objects in the scene.
// If the ray misses all objects, it returns black (`Vec()`).
// If it hits an object, it computes the radiance based on the material type (diffuse, specular, or refractive).

// `Ray r`**: The incoming ray that is being traced.
// `int depth`: The recursion depth. This controls how many times a ray can reflect or refract before the path tracing stops. 
// If `depth` exceeds a certain value (in this case, 5), the algorithm terminates early with a probability.
// `unsigned short *Xi`**: A pointer to a random number generator state array. 
// It's used for randomness (sampling rays) to simulate the stochastic nature of light transport.
Vec radiance(const Ray &r, int depth, unsigned short *Xi){


  // The first step is checking if the ray `r` intersects with any objects in the scene (which are stored in the `spheres[]` array). 
  // The `intersect` function returns the closest intersection point (`t`) and the ID of the intersected object (`id`).
  // If there is no intersection (i.e., the ray misses everything), the function returns `Vec()`, which represents black (no light).
  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return Vec(); // if miss, return black

  const Sphere &obj = spheres[id];        // The object hit by the ray
  // The `obj` is the sphere that the ray intersected. It is used to access the color (`c`), emission (`e`), and reflection type (`refl`) of the object, as well as to calculate the further behavior of the ray.

  // `x` is the point of intersection on the sphere.
  // `n` is the normal vector at that point, which is the vector from the sphere's center to the intersection point, normalized.
  // `nl` is the adjusted normal direction. If the ray is coming from inside the sphere, the normal should point inward (multiply by -1).
  Vec x = r.o + r.d * t;         // Intersection point
  Vec n = (x - obj.p).norm();    // Normal at the intersection point
  Vec nl = n.dot(r.d) < 0 ? n : n * -1; // Adjust normal direction
  Vec f=obj.c;

  double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl

   // This block controls recursion depth using Russian Roulette (RR). 
   // If the depth exceeds 5, there's a chance (based on the surface's reflectivity) that the function will terminate early.
   // If the algorithm chooses to stop, the function will return the sphere's emission color (`obj.e`), simulating light sources in the scene.
  if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return obj.e; //R.R.
  

  // In this case, the surface has a diffuse reflection (e.g. a matte surface).
  // The function generates a random direction `d` in 3D space that is scattered around the normal `n`.
  // It then calls the `radiance` function recursively with the new ray (`Ray(x, d)`), and the final color is the emission color `obj.e` plus the reflected color (`f.mult(...)`).
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
    Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
    return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
  // In the case of **specular reflection** (e.g. a mirror), the ray bounces off the surface in a perfect reflection, 
  // so the new direction is calculated as the reflection of the incoming ray (`r.d - n * 2 * n.dot(r.d)`).
  // The function then recursively calls `radiance` to trace the reflected ray. 
  } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));

  // If the material is **refractive** (e.g. glass), the function calculates the refraction direction `tdir` based on Snell's law 
  // and handles both reflection and refraction with **total internal reflection** (TIR).
  // Similar to diffuse and specular, the function recursively traces the refracted ray (`Ray(x, tdir)`).
  // The function uses **Russian Roulette** (RR) to probabilistically choose between reflection 
  // and refraction based on the reflection coefficient (`Re`), giving a more physically accurate result.
  
  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION

  bool into = n.dot(nl)>0;                // Ray from outside going in?

  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;

  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + f.mult(radiance(reflRay,depth,Xi));

  Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();

  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));

  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);

  return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

// The `main` function sets up the camera and image dimensions, performs supersampling to reduce aliasing, 
// and traces rays through the scene to compute radiance for each pixel using recursive path tracing. 
// Command-line arguments** (`argc` and `argv`) are used to determine how many samples per pixel should be used for supersampling.
int main(int argc, char *argv[]){

  // `w` and `h`: These define the width (`1024` pixels) and height (`768` pixels) of the image.

  // `samps`: This defines the number of samples per pixel for supersampling. 
  // If the user provides a second command-line argument (`argc == 2`), it is used to determine the number of samples. 
  // This value is divided by 4, which seems to be a performance optimization. Otherwise, the default is `1` sample.

  // `cam`: This defines the camera’s position and direction. 
  // The camera is located at `(50, 52, 295.6)` and looks toward the direction `(0, -0.042612, -1)`. 
  // The `norm()` function normalizes the direction vector so it has a magnitude of 1.

  // `cx` and `cy`: These vectors define the horizontal and vertical direction vectors of the camera’s view plane. 
  // They are scaled based on the image resolution (`w`, `h`), and `cx` is adjusted based on the aspect ratio. 
  // `cy` is computed as the perpendicular vector to `cx` in the camera's coordinate system. 
  // These vectors help to map screen space coordinates to 3D world space.
  
  // `c`: This is an array of `Vec` objects that will store the RGB color of each pixel in the image. 
  // The size of this array is `w * h`, and each pixel is initialised to zero.

  int w=1024, h=768, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples

  Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir

  Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];

// `#pragma omp parallel for`: This line enables parallel processing using OpenMP. 
// The `for` loop that follows is parallelized, meaning each row (`y`) of pixels is computed in parallel by multiple threads, improving performance.
// `schedule(dynamic, 1)` means that OpenMP will distribute the rows of pixels dynamically among threads. Each thread will process one row at a time and may pick the next available row to process as soon as it finishes a previous one.
// `private(r)` ensures that each thread gets its own local variable `r` (used to accumulate radiance) during rendering.
  
  #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y=0; y<h; y++){                       // Loop over image rows
      fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
      for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols
        for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
          for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
            
// Outer loop (`y`): Iterates over each row of pixels (`y` axis).
// Inner loop (`x`): Iterates over each column of pixels (`x` axis).
// Sub-pixel sampling (`sy`, `sx`): This part implements 2x2 supersampling for each pixel to reduce aliasing.
// The pixel area is divided into a 2x2 grid, and for each subpixel (sy, sx), the radiance is computed.
            for (int s=0; s<samps; s++){

// `r1` and `r2`: These are random numbers generated using `erand48(Xi)`. 
// These values are used to add a bit of randomness to the subpixel sampling (for anti-aliasing purposes).
// `dx` and `dy`: These values are computed based on `r1` and `r2` to generate jittered offsets for the subpixel sampling within each pixel. They are used to adjust the ray's direction slightly so that multiple rays sample different parts of the pixel.

// `Vec d`: This is the direction vector for the ray. 
// The formula calculates the direction from the camera to the pixel, adding in the jitter from `dx` and `dy`. 
// `cx` and `cy` scale this direction to the correct position on the image plane. 
// Then, the camera's direction `cam.d` is added to get the final ray direction.

// `radiance(Ray(cam.o+d*140, d.norm()), 0, Xi)`: This is the key function that calculates the radiance (or color) of the pixel. 
// It generates a ray that starts at the camera’s position (`cam.o`) and moves in the direction `d`. 
// The `140` factor scales the ray to simulate a large view distance, making sure that the rays travel far enough into the scene. 
// The `radiance` function is called recursively, calculating the radiance along the path of the ray through the scene.

// The value returned by `radiance` is added to the accumulator `r`, 
// and the result is averaged by dividing by `samps` to account for the number of samples per pixel.

              double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
              double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
              Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                      cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
              r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
            } // Camera rays are pushed ^^^^^ forward to start in interior
// `clamp(r.x), clamp(r.y), clamp(r.z)`: This ensures that the color values are within the range `[0, 1]`, which is important for color normalization.
// `*.25`: This is a simple gamma correction (approximated by multiplying by `.25`). 
// It's a form of tone mapping that brightens the colors for display, 
// as the radiance values can be much larger than 1 due to the accumulation of light from multiple samples.
            c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
          }
  }

// Opening the file and saving the image
  FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
  fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
