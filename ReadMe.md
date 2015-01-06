Ray Trace Sample
==========================================

This program is a simple Ray Tracer. It is written in C/C++ with GLUT. While it is mostly completed,
There are a few things that are missing. The first is that there is a bug somewhere in displaying
the Bspline Triangle Mesh. This mesh (which has it's own separate program entitled BsplineMeshMaker.cpp)
has a problem in it's program in that it does not correctly display the texture given. This can be, however,
entirely ignored in the ray tracer as all the ray tracer requires from the mesh are the points which define it.
The second property that is missing is that it completely does not compute the any refraction ray.

I may or may not return to work out the bugs as I am currently interested in other projects at the moment,
however, this does show many concepts needed to inspire a working Ray Tracer written in C/C++ with OpenGl.
This was written in Visual Studio 2013.
