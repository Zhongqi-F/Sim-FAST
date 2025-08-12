clear all;
clc;
close all;

tic

%% Define the geometry of experiment samples
L=305*10^(-3)/2;
W=80*10^(-3);
alpha=180*pi;
t=20*10^(-3); 
gap=0*10^(-3); 


% Stiffness parameters of the structure
faceThickness=2*10^(-3); 
E=4*10^9; % Young's modulus
v=0.3; % Posson's Ratio


% Initialize Elements
node=Elements_Nodes; 
rot_spr_4N=Vec_Elements_RotSprings_4N_Directional; 
cst=Vec_Elements_CST; 
zlspr=Vec_Elements_Zero_L_Spring; 

assembly=Assembly_ThickOrigami();
assembly.node=node;
assembly.cst=cst;
assembly.rot_spr_4N=rot_spr_4N;
assembly.zlspr=zlspr;


%% Define the nodal coordinates
node.coordinates_mat=[0 0 0;
                      2*L/3 0 0; %% NODE 2 
                      L 0 0; % NODE 3
                      L W 0;
                      2*L/3 W 0; %% NODE 5
                      0 W 0; % NODE 6
                      0 0 t;
                      2*L/3 0 t; %% NODE 8
                      L 0 t; % NODE 9
                      L W t;
                      2*L/3 W t; %% NODE 11
                      0 W t; % NODE 12
                      gap+L 0 0;
                      gap+4*L/3 0 0; %% NODE 14
                      gap+2*L 0 0; %NODE 15
                      gap+2*L W 0;
                      gap+4*L/3 W 0; %% NODE 17
                      gap+L W 0; % NODE 18
                      gap+L 0 t;
                      gap+4*L/3 0 t; %% NODE 20
                      gap+2*L 0 t; % NODE 21
                      gap+2*L W t;
                      gap+4*L/3 W t; %% NODE 23
                      gap+L W t];  % NODE 24
                      

% set up panels for plotting
% sample 1
cst.node_ijk_mat=[1 2 8;  %% cst.node_ijk_mat
                 1 7 8;
                 2 3 8;
                 3 8 9;
                 3 4 10;
                 3 9 10;
                 4 10 11;
                 4 5 11;
                 5 6 11;
                 6 11 12;
                 1 6 12;
                 1 7 12;
                 7 8 12;
                 8 11 12;
                 8 9 10;
                 8 11 10;
                 1 2 5;
                 1 5 6;
                 2 3 5;
                 3 4 5];

% sample 2
cst.node_ijk_mat=[cst.node_ijk_mat; % this is for double set up
                 13 14 20;
                 13 19 20;
                 14 15 20;
                 15 20 21;
                 15 16 22;
                 15 21 22;
                 16 17 22;
                 17 22 23;
                 17 18 24;
                 17 23 24;
                 13 18 24;
                 13 19 24;
                 19 20 23;
                 19 23 24;
                 20 21 23;
                 21 22 23;
                 13 14 17;
                 13 17 18;
                 14 15 17;
                 15 16 17];


rot_spr_4N.rot_spr_K_vec=1; 
rot_spr_4N.mv_factor_vec=1000;

cst.t_vec=faceThickness*ones(40,1); %%%%%
cst.v_vec=v*ones(40,1);
cst.E_vec=E*ones(40,1);

% Springs for the panel bending (locked)
zlspr.node_ij_mat=[10 24;
                   9 19;
                   4 18; 
                   3 13];

% hinge down
% rotSpr.node_ijkl_mat=[2 3 4 17]; % Set up the rotational hinges
% rotSpr.mv_vec=0;
% zlspr.k_vec=[119999;
%     119999; 
%     215806;  
%     215806]; 
% load=-100/50;


% hinge up 
rot_spr_4N.node_ijkl_mat=[8 9 10 23]; % Set up the rotational hinges
rot_spr_4N.mv_vec=1;
zlspr.k_vec=[215806; 
    215806; 
    119999;  
    119999]; 
load=-40/50;

assembly.Initialize_Assembly();

plots=Plot_ThickOrigami;
plots.assembly=assembly;
plots.displayRange=[-0.1;0.4;-0.1;0.2;-0.1;0.2];
plots.Plot_Shape_NodeNumber;
plots.Plot_Shape_ZLsprNumber;

%% Set up the loading solver
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.supp=[1,1,1,1; % NODE 1
         6,1,1,1; % NODE 6
         15,0,1,1;
         16,0,1,1;]; 


nr.tol=10^-7;
nr.increStep=50;

% Set up the load (applied at 4 nodes)
nr.load=[8,0,0,load;
         11,0,0,load;
         20,0,0,load;
         23,0,0,load;];

% Solve for the deformation history
Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));

%% Plotting results
nodes_of_interest = [8, 11, 20, 23];
num_steps = size(Uhis, 1);

displacements_mm = zeros(num_steps, length(nodes_of_interest));
for i = 1:length(nodes_of_interest)
    displacements_mm(:, i) = abs(Uhis(:, nodes_of_interest(i), 3)) * 1000;
end
displacements_mm_average=mean(displacements_mm,2);
displacements_mm_average=[0;displacements_mm_average];

reaction_load_N = linspace(0, 4*abs(load)*50, num_steps+1)'; 
reaction_load_kN = reaction_load_N / 1000; 

figure;
hold on;
plot(displacements_mm_average, reaction_load_kN);
hold off;

xlabel('Displacement (mm)');
ylabel('Load (kN)');
title('Displacement vs Load (hinge up)');
% legend('Location','best');
grid on;
box on;

toc