3D Graffiti
==========================

This OpenGL demo app reads BVH motion capture data taken from a hand held 
controller and generates a triangle mesh from the path. Controller rotation data
is read but is currently ignored as I found it didn't give results as good as the
parallel transport method for constructing a reference frame along a path.
A nice description of the paralle transport method is given in 
"Quaternion Gauss Maps and Optimal Framings of Curves and Surfaces (1998)", 
http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.50.8658).

Lighting is performed using a third order spherical harmonic ambient
term with a direct spot light and PCF filtered shadow mapping. The ambient
SH coefficients are calculated from Paul Debevec's HDR "beach" light probe.

The mesh is reconstructed from scratch each frame and nothing is
optimised so I expect it could be made significantly faster.

Controls
--------

mouse - Rotate camera
w,a,s,d - Translate camera

1 - Open Animation (.bvh)
2 - Open Mouse Track (.txt)
3 - Export mesh at current frame (creates dew.obj)
4 - Start/Stop Mouse Record (saves mouse clicks to control.txt)

+/- Adjust Mouse track sync forward/back in time

j - Speed up animation
u - Slow down animation
r - Rewind to beginning
f - Fast forward to end
space - Pause animation

i - Increase path width
k - Decrease path width
o - Increase path thickness
l - Decrease path thickness
p - Increase the effect velocity has on path thickness
; - Decrease the effect velocity has on path thickness
/ - Cycle brushes

h - Hide help text
c - Hide can
b - Show wirefame 
n - Show normals
b - Toggle fullscreen
