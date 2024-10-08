## Simulator For Active STructures (Sim-FAST)

This is an educational and research code package for active structures, including tensegrity, 
origami, trusses, mechanisms, metamaterials, and many others. **Please ensure running in MATLAB 2024a or newer version!**

The following Intro Figure shows some working examples.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/Figures/Intro.png)

More specifically, the model in this package represents an active structure 
using nodal representations without rotational degrees-of-freedom. 
Such formulation are similar to the pseudo-rigid-body model coined in the 
pioneering research by Prof. Larry Howell. 
The ability to capture phenonmena like bending and twisting is enabled 
through using four-node and three-node rotational spring elements. 

While the modeling approach may not be suitable for all active structures 
(there probably exists no universal model for every problem), I believe the 
framework is particularly suitable for educational purpose and many 
interesting research problems I myself have worked on. 
This formulation is good for demonstrating element formulation for
large deformation, implicit and explicit nonlinear solvers, 
and integration of simulation package. 
The formulation can be easily adopted for many active structures rapidly
for a research project, giving people great flexibility. 
Being non-black-box also make it easy to make changes to the code for 
specific problems or connect the code with optimization packages for design. 
I will not claim that the code package is a replacement to commercial 
softwares (which is not the case), but it offers very different capabilities 
that I myself does not find in commercial softwares. 

 

## A Living Textbook 

This package is associated with  a living textbook component. 
This living textbook is based on the live script function in MATLAB, where we can combine text, equations, figures, codes, functions, and results in one place. 
This is really ideal for introducing simulation of active structures because we really need all these components side by side. 
I am constructing this package for my own students, I hope this will
become their starting point to learn simulation for active structures. 
However, I believe such an effort is worth made open-access, which is 
why I uploaded this package here.

![alt text](https://github.com/zzhuyii/Sim-FAST/blob/main/Figures/Feature.png)

**Please ensure running in MATLAB 2024a or newer version!**

In the live script ("living textbook"), functions are created as we are working through examles. Thus, deifinition of functions are merged with global codes. Running such mixture of live script is only supported in 2024a or newer. In older MATLAB versions, functions must be in the back (which is not ideal when the code is integrated into textbooks). This "living-textbook" idea is only possible with the 2024a or newer MATLAB so Hey! We are definitely the state-of-the-art!

When using the code, please add "00_SourceCode_Elements" and "00_SourceCode_Solver" 
to the path in MATLAB. For some examples, we can use "00_SourceCode_Elements_Advanced" 
for faster computation. 

Basic elements are derived using central-difference method while 
advanced elements are derived analytically. Advanced elements are all vectorized for
computation speed and they begin with "Vec_Elements". Basic elements begin with "CD_Elements"
representing that their derivation is based on central difference. 





