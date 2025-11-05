# ray-tracing-scene-builder
A simple scene builder sub-module for ray tracing module  

Using Smallpt code as a base (https://www.kevinbeason.com/smallpt/).

* smallpt.cpp - Original code with recursion
Time: 1m29.429s
* explicit.cpp - Speedup for small lights. Adds explicit light sampling and a small function signature change. No speecup for large lights like in the default image. 
Time: 2m1.625s
* forward.cpp - Revision of radiance() function that removes all recursion and uses only a simple loop and no path branching. That is, the ray tree is always one ray wide. Does not provide speedup and in fact makes it slower. This implementation has the same disadvantages as the original code with recursion, as it uses while (1) with not parameter explicitly stopping loops.
Time: 3m23.045s
