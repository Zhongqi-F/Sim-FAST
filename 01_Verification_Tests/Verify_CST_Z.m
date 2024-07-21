%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% Here we generate the geometry of the origami

% Define the nodal coordinate before meshing
L=100*10^(-3);
W=100*10^(-3);

% Stiffness parameters of the structure
sprStiff=0.01;

panelE=10*10^6;
panelv=0.3;
panelt=1*10^-3;

node=Elements_Nodes;

cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings;



%% Geometry of a Plate
% nodal coordinates
node.coordinates_Mat=[
    0,0,0;
    0,L,0;
    W,0,0;
    W,L,0];


% set up panels for plotting
cst.cst_ijk_mat=[1,2,4;
                 1,3,4];

cst.v_vec=panelv*ones(2,1);
cst.E_vec=panelE*ones(2,1);
cst.t_vec=panelt*ones(2,1);


%% Set up the rotational springs
rotSpr.rotSprIJKL_Mat=[3 1 4 2];

rotSpr.rotSprK_Vec=sprStiff;


%% Initialize assembly
assembly=Assembly_CST_Origami();
assembly.node=node;
assembly.cst=cst;
assembly.rotSpr=rotSpr;

assembly.InitializeAssembly()


%% Test CST

delta=10^-4;
U=[0,0,0;
   0,delta,0;
   0,0,0;
   0,delta,0;];

[F,K]=cst.SolveFK(node,U);


% [bar_strain_mat,l_mat,x_reshape]=...
%     cst.Calc_Bar_Strain(U,node.coordinates_Mat);
% 
% [cst_strain_mat] = cst.Calc_CST_Strain(bar_strain_mat);
% 
% [dedx,d2edx2] = cst.Calc_Derivatives(x_reshape,l_mat);
% 
% [Tcst]=cst.GlobalForce(U,dedx,cst_strain_mat);
% 
% [Kcst]=cst.GlobalStiff(U,dedx,d2edx2,cst_strain_mat);


%% Plot for investigation
plots=Plot_CST_Origami();
plots.displayRange=0.15;
plots.displayRangeRatio=0.4;
plots.assembly=assembly;

plots.Plot_Shape_CSTNumber()
plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_SprNumber()



%% Setup the loading controller
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.supp=[1,1,1,1;
         3,1,1,1;
         4,1,1,1;];

force=0.01;
nr.load=[2,0,0,force];

nr.increStep=50;
nr.tol=10^-5;
nr.iterMax=50;

Uhis=nr.Solve();

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));

% plots.fileName='CST_Kresling.gif';
% plots.Plot_DeformedHis(Uhis);

%% Find the reaction force and loading results

forceHis=zeros(nr.increStep,1);
UrefHis=zeros(nr.increStep,1);

for i=1:nr.increStep
    [F,K]=assembly.SolveFK(squeeze(Uhis(i,:,:)));
    UrefHis(i)=Uhis(i,2,3);
    forceHis(i)=F(2*3);
end

figure
UrefHis=[0;UrefHis];
forceHis=[0;forceHis];
plot(UrefHis,forceHis)
xlabel('Deformation (m)') 
ylabel('Applied Force (N)') 

