## Simulator For Active STructures (Sim-FAST)

This is an educational code package for active structures, including tensegrity, 
origami, trusses, mechanisms, metamaterials, and many others. 
As of now, this package is still under construction. 
Eventually, there will be a note associated with this package. 
Hoepfully I will get the opportunity of teaching it one day. 

The following Intro Figure shows a few working examples.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/Figures/Intro.png)

More specifically, the model in this package represents an active structure 
using nodal representations without rotational degrees-of-freedom. 
Such formulation are similar to the pseudo-rigid-body model coined in the 
pioneering research of Prof. Larry Howell. 
The ability to capture phenonmena like bending and twisting is enabled 
through using four-node and three-node rotational spring elements. 

While the modeling approach may not be suitable for all active structures 
(there probably exists no universal model for every problem), I believe the 
framework is particularly suitable for showing element formulation for large
deformation, implicit and explicit nonlinear solvers, and integration of 
simulation package to a student trying to study active structure. 
This is why I started this package here. 

## As of Jan 13th. 

The code is working for a few examples, including active truss, thin origami, 
thick origami, mechanism, and knitting structures. 
Not particularly sure how much effort I can spend on this as the semester start. 
Will continue working on the code when have time. 


## Using the code

When using the code, please add "00_SourceCode_Elements" and "00_SourceCode_Solver" 
to the path in MATLAB. For some examples, we can use "00_SourceCode_Elements_Advanced" 
for faster computation. However, advanced elements (vectorized) are not available for
every example for now. 





