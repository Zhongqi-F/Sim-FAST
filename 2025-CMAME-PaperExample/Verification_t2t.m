clear all
close all
clc

tic

%% Define Geometry
node=Elements_Nodes;
randomNodeMat=rand([3 2]);

node.coordinates_mat=[node.coordinates_mat;
    0 0 0;
    0 1 0;
    1 0 0;
    randomNodeMat  0.5*ones(3,1)+0.4*rand([3,1]);
    randomNodeMat  2*ones(3,1);];


%% Define assembly
assembly=Assembly_MembraneWithBar;
cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
t2t=CD_Elements_T2T_Contact;
bar=CD_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.t2t=t2t;
assembly.rotSpr=rotSpr;
assembly.bar=bar;

%% Define Plotting Functions
plots=Plot_MembraneWithBar;
plots.assembly=assembly;
plots.displayRange=[-0.5; 1.5; -0.5; 1.5; -0.5; 1];

plots.viewAngle2=20;

plots.Plot_Shape_NodeNumber;


%% Define Triangle
cst.node_ijk_mat=[
    1 2 3;
    4 5 6;];

% material properties
t=0.001;
E=10^9;
cst.t_vec=t*ones(2,1); % Thickness 
cst.E_vec=E*ones(2,1); % Young's Modulus
cst.v_vec=0.2*ones(2,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;
plots.Plot_Shape_NodeNumber;

%% Define Rotational Spring
rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
    1 2 3 4];

rotSpr.rot_spr_K_vec=0;
rotSpr.theta_stress_free_vec=pi;
% Create a fake spring

 
%% Define Triangle to Triangle Penetration Prevention
t2t.d0=0.02;
t2t.delta=10^-5;
t2t.k_contact=1;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[1;2];

% To deactivate penetration prevention
% t2t.group_number=[1;1];


%% Set up bar
bar.node_ij_mat=[4 7;
                 5 8;
                 6 9;];
bar.E_vec=ones(3,1);
bar.A_vec=0.1*ones(3,1);

% set up ploting color
plots.colorNum=[1;2;];

assembly.Initialize_Assembly;



%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[1 1 1 1;
     2 1 1 1;
     3 1 1 1;
     7 1 1 1;
     8 1 1 1;
     9 1 1 1;
     4 1 1 0;
     5 1 1 0;
     6 1 1 0;];

force=0.004;
step=50;

Uhis=[];
for k=1:step
    nr.load=[6 0 0 -force*k;];
    
    nr.increStep=1;
    nr.iterMax=5;
    nr.tol=1*10^-6;

    Uhis(k,:,:)=squeeze(nr.Solve());
end

toc
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='Verification_t2t.gif';
plots.Plot_DeformedHis(Uhis)

