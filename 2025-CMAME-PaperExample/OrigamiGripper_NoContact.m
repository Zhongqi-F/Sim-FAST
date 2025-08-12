clear all
close all
clc

tic

%% Define Geometry
L1=0.001;
L2=0.0003;
L3=0.0004;
L4=0.0005;
La=0.0007;

W=0.0003;
Wh=0.00012;
Wt=0.00004;

node=Elements_Nodes;

node.coordinates_mat=[node.coordinates_mat;
    0   0    0;
    -W   W/4  0;
    W    W/4  0;
    -W   W/4+L3  0;
    W    W/4+L3  0;  %5
    -W   W/4+L3+L4  0;
    W    W/4+L3+L4  0;
    -Wh  W/4+L3  0;
    Wh    W/4+L3  0;
    -Wh  W/4+L3+L4/3  0;   %10
    Wh    W/4+L3+L4/3  0;
    -W-La  W/4+L3  0;
    -W-La  W/4+L3+L4/2  0;
    W+La  W/4+L3  0;
    W+La  W/4+L3+L4/2  0; %15

    0   W/2+L3+L4    0;
    -W  3*W/4+L3+L4    0;
    W   3*W/4+L3+L4    0;
    -W  3*W/4+L3+L4+L1    0;
    -Wt   3*W/4+L3+L4+L1    0; %20
    Wt  3*W/4+L3+L4+L1    0;
    W   3*W/4+L3+L4+L1    0;
    -Wt   3*W/4+L3+L4+L1+L2    0; 
    Wt  3*W/4+L3+L4+L1+L2    0;
    ];




%% Define assembly
assembly=Assembly_Membrane;
cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
t2t=CD_Elements_T2T_Contact;

assembly.cst=cst;
assembly.node=node;
assembly.t2t=t2t;
assembly.rotSpr=rotSpr;


%% Define Plotting Functions
plots=Plot_Membrane;
plots.assembly=assembly;
plots.displayRange=[-0.002; 0.002; -0.001; 0.004; -0.001; 0.003];

plots.Plot_Shape_NodeNumber;


%% Define Triangle

cst.node_ijk_mat=[cst.node_ijk_mat;
    1 2 3;
    2 4 8;
    2 8 9;
    2 3 9;
    3 5 9;
    12 4 13;
    4 6 13;
    4 6 10;
    4 8 10;
    6 7 10;
    10 11 7;
    9 11 5;
    5 7 11;
    5 14 15;
    5 7 15;
    16 17 18;
    17 19 20;
    17 20 21;
    17 18 21;
    21 22 18;
    20 21 23;
    21 23 24;];
 
% thickness of the triangle
t=0.00002;
E=2*10^9;

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);

cst.t_vec=t*ones(cstNum,1); % 1mm thickness
cst.E_vec=E*ones(cstNum,1); % 10^8 Young's Modulus
cst.v_vec=0.2*ones(cstNum,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;

%% Define Rotational Spring

rotSpr.node_ijkl_mat=[
    1 2 3 9;
    3 2 9 8;
    9 2 8 4;
    2 4 8 10;
    8 4 10 6;
    10 4 6 13;
    12 4 13 6;
    3 5 9 11;
    9 5 11 7;
    11 5 7 15;
    7 5 15 14;
    5 11 7 10;
    11 7 10 6;
    7 6 10 4;
    16 17 18 21;
    18 17 21 20;
    21 17 20 19;
    17 20 21 23;
    20 21 23 24;
    17 21 18 22;
    ];


% Find the bending stiffness
rotk=E*(t*3)^3/6/sqrt(3)/(1-0.2^2);
rotNum=size(rotSpr.node_ijkl_mat);
rotNum=rotNum(1);

for i=1:rotNum
    node1=rotSpr.node_ijkl_mat(i,2);
    node2=rotSpr.node_ijkl_mat(i,3);
    Lbar=norm(node.coordinates_mat(node1,:)-node.coordinates_mat(node2,:));
    rotSpr.rot_spr_K_vec=rotk*ones(rotNum,1);
end

% we assume that the panels are 3 time thicker than the folding line
rotSpr.rot_spr_K_vec([1,6,10,15])=rotSpr.rot_spr_K_vec([1,6,10,15])/27;
rotSpr.theta_stress_free_vec=pi*ones(rotNum,1);

plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
t2t.d0=0.00004;
t2t.delta=10^-9;
t2t.k_contact=100;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[ones(15,1);
                  ones(7,1);];
% Set group to be all 1 to avoid contact activation


% Set up color for plotting
plots.colorNum=[ones(15,1);
                2*ones(7,1);];
assembly.Initialize_Assembly;

%% Set up solver
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;

sf.supp=[
    1 1 1 1;    
    2 1 1 1;
    3 1 1 1;
    16 1 1 1;
    17 1 1 1;
    18 1 1 1;
];

sf.increStep=80;
sf.tol=10^-4;
sf.iterMax=10;

sf.increStep=20;
sf.targetRot=rotSpr.theta_stress_free_vec;
sf.targetRot(1)=pi+0.6*pi;
Uhis1=sf.Solve;

sf.increStep=20;
sf.targetRot(15)=pi+0.89*pi;
Uhis2=sf.Solve;

sf.increStep=80;
sf.targetRot(1)=pi+0*pi;
Uhis3=sf.Solve;

sf.increStep=80;
sf.targetRot(6)=pi+0.5*pi;
sf.targetRot(10)=pi+0.5*pi;
Uhis4=sf.Solve;

toc
plots.Plot_DeformedShape(squeeze(Uhis1(1,:,:)))
plots.Plot_DeformedShape(squeeze(Uhis1(end,:,:)))
plots.Plot_DeformedShape(squeeze(Uhis2(end,:,:)))
plots.Plot_DeformedShape(squeeze(Uhis3(end,:,:)))
plots.Plot_DeformedShape(squeeze(Uhis4(end,:,:)))
 
Uhis=cat(1,Uhis1,Uhis2);
Uhis=cat(1,Uhis,Uhis3);
Uhis=cat(1,Uhis,Uhis4);

plots.fileName='OrigamiGripper_NoContact.gif';
plots.Plot_DeformedHis(Uhis)

