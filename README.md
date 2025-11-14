# ray-tracing-scene-builder
A simple scene builder sub-module for ray tracing module  

Using Smallpt code as a base (https://www.kevinbeason.com/smallpt/).

1. smallpt.cpp - Original code with recursion
Time (g++-13 -std=c++20 -Ofast -march=native -fopenmp main.cpp -o main): 
0m2.718s
2. explicit.cpp - Speedup for small lights. Adds explicit light sampling and a small function signature change. No speecup for large lights like in the default image. 
Time (g++-13 -std=c++20 -Ofast -march=native -fopenmp main.cpp -o main):  
0m3.672s
3. forward.cpp - Revision of radiance() function that removes all recursion and uses only a simple loop and no path branching. That is, the ray tree is always one ray wide. Does not provide speedup and in fact makes it slower. This implementation has the same disadvantages as the original code with recursion, as it uses while (1) with not parameter explicitly stopping loops.
Time (g++-13 -std=c++20 -Ofast -march=native -fopenmp main.cpp -o main):  
0m3.980s
4. main1.cpp - Substituted vec structure with Eigen::Vector3d from Eigen C++ library.
Time (g++-13 -std=c++20 -Ofast -march=native -fopenmp main.cpp -o main):
0m2.817s

main2.cpp - Added the option to create a triangle in the scene
main3.cpp - Added .obj file processing


Compilation:
g++-13 -std=c++20 -Ofast -march=native -fopenmp main.cpp -o main
-march=native: Optimise code for the host machine's specific architecture
-fopenmp: Enable OpenMP support
-DNDEBUG: Disable debugging-related code or assertions during the compilation
std=c++20
Optimisation flags:
-O1, O2, O3: Various levels of optimisation
-Ofast: It enables a range of optimisation techniques that improve the performance of the code. This typically includes all optimisations performed by -O2, along with additional optimisations that can improve runtime speed but may break certain assumptions of the language standard or introduce undefined behavior in some cases.
Problem1: O1, O2, O3, Ofast produce completely incorrect outputs. 
Solution1: r vector was not initialised to (0, 0, 0) at the beginning of each loop.



Original smallpt.cpp modifications:
*  Corrected this warning
"main.cpp: In function ‘int main(int, char**)’:
main.cpp:123:46: warning: narrowing conversion of ‘((y * y) * y)’ from ‘int’ to ‘short unsigned int’ [-Wnarrowing]
  123 |       for (unsigned short x=0, Xi[3]={0,0,y*y*y}; x<w; x++)   // Loop cols"
* Got rid of the vector structure and replace it with the standard library (e.g. Eigen)
* Changed some if-else structures from one liners to a readable indented structure
* Added the option to create a triangle in the scene
* Added .obj files import and processing
* Amended camera settings to add Field of View (FOV)
* Refactored the code to add intersectable class
* Refactored the code to add intersectable class


To do:
* Remove recursion (try stacks)
* Implement Bounding Volume Hierarchy (BVH)
* Camera settings
* Check how to verify the output correctness
* Try it on deformation example

"Ray Casting Curved-Quadratic Elements" paper

Some links:
https://libeigen.gitlab.io/eigen/docs-nightly/GettingStarted.html
https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html


