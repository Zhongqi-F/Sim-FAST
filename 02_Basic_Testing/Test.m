clear all
close all
clc
tic

%% Define Geometry
node=Elements_Nodes;
node.coordinates_mat=[0 0 0;
                      0 1 1;
                      0 0 2;
                      0 -1 1;
                      0 0 1;];

disp1=0;
disp2=0.5;

nodispMat=[0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;
           0 0 0;];

dispMat=[0 0 0;
         disp1 disp2 0;
         0 0 0;
         0 0 0;
         0 0 0;];

%% Define assembly
assembly=Assembly_Kirigami_Truss;
cst=Vec_Elements_CST;
rot_spr_4N=Std_Elements_RotSprings_4N;
rot_spr_3N=Std_Elements_RotSprings_3N;
bar=Std_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_4N=rot_spr_4N;

%% Define Plotting Functions
plots=Plot_Kirigami_Truss;
plots.assembly=assembly;
plots.displayRange=[-2; 2; -2; 2; -2; 2];

plots.Plot_Shape_Node_Number;


%% Define Bars
bar.node_ij_mat=[1 2;
                 2 3;];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=ones(barNum,1);
bar.E_vec=ones(barNum,1);

plots.Plot_Shape_Bar_Number();

%% Define Triangle
cst.node_ijk_mat=[1 2 3;
                  1 3 4];

cstNum=size(cst.node_ijk_mat,1);
cst.t_vec=ones(cstNum,1);
cst.E_vec=ones(cstNum,1);
cst.v_vec=0.2*ones(cstNum,1);

plots.Plot_Shape_CST_Number;



%% Define Rotational Spring
rot_spr_4N.node_ijkl_mat=[1 2 4 3;
                        2 1 3 4;];

rotNum=size(rot_spr_4N.node_ijkl_mat);
rotNum=rotNum(1);
rot_spr_4N.rot_spr_K_vec=ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;

assembly.Initialize_Assembly;

%% Define 3 Node spring

rot_spr_3N.node_ijk_mat=[1 2 3;
                         1 4 3;
                         1 5 3];
rot_spr_3N.rot_spr_K_vec=[1;1;1];

rot_spr_3N.Initialize(node);



%% Get forces and stiffness

[F,K]=assembly.Solve_FK(dispMat);
[Fcst,Kcst]=cst.Solve_FK(node,dispMat);
[Fbar,Kbar]=bar.Solve_FK(node,dispMat);
[Frs4,Krs4]=rot_spr_4N.Solve_FK(node,dispMat);
[Frs3,Krs3]=rot_spr_3N.Solve_FK(node,dispMat);
