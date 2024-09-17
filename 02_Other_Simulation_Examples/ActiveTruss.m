clear all;
clc;
close all;

%% Define the geometry of the active truss
% Define the nodes
t=0.5;
node=Elements_Nodes;
node.coordinates_mat=[0 0 t; 
                      1 0 t;
                      2 0 t; 
                      3 0 t; 
                      4 0 t; % 5
                      0 1 t; 
                      1 1 t;
                      2 1 t; 
                      3 1 t; 
                      4 1 t; % 10
                      0.5 0.5 2*t; 
                      1.5 0.5 2*t;
                      2.5 0.5 2*t; 
                      3.5 0.5 2*t; 
                      0.5 0.5 0; 
                      1.5 0.5 0;
                      2.5 0.5 0; 
                      3.5 0.5 0; ];

% Define the bars
bar=CD_Elements_Bars;

% The bars for the square base
bar.node_ij_mat=[1,2;
                    2,3;
                    3,4;
                    4,5;
                    6,7;
                    7,8;
                    8,9;
                    9,10;
                    1,6;
                    2,7;
                    3,8;
                    4,9;
                    5,10];

% Top and bottom row
bar.node_ij_mat=[bar.node_ij_mat;                    
                    11,12;
                    12,13;
                    13,14;
                    15,16;
                    16,17;
                    17,18;];

% Diagonal bars
for i=[1,2,6,7]
    bar.node_ij_mat=[bar.node_ij_mat;11,i];
    bar.node_ij_mat=[bar.node_ij_mat;15,i];
end

for i=[2,3,7,8]
    bar.node_ij_mat=[bar.node_ij_mat;12,i];
    bar.node_ij_mat=[bar.node_ij_mat;16,i];
end

for i=[3,4,8,9]
    bar.node_ij_mat=[bar.node_ij_mat;13,i];
    bar.node_ij_mat=[bar.node_ij_mat;17,i];
end

for i=[4,5,9,10]
    bar.node_ij_mat=[bar.node_ij_mat;14,i];
    bar.node_ij_mat=[bar.node_ij_mat;18,i];
end
    
%% Initialize area and other properties
bar.A_vec=ones(51,1);
bar.E_vec=ones(51,1);

bar.Initialize(node);

% Create Assembly
assembly=Assembly_Truss;
assembly.bar=bar;
assembly.node=node;
assembly.InitializeAssembly();


plots=Plot_Truss();
plots.displayRangeRatio=0.1;
plots.displayRange=5;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()

%% Set up the loading solver
action=Solver_NR_TrussAction;
action.assembly=assembly;
action.supp=[1,1,1,1;
         6,1,1,1;
         15,1,1,1;];

% Set up the stress-free length of the bars
action.targetL0=bar.L0_vec;
dL0=0.2;

% Set up how the stress-free lengths of bars change
action.targetL0(14)=action.targetL0(14)-dL0;
action.targetL0(15)=action.targetL0(15)-dL0;
action.targetL0(16)=action.targetL0(16)-dL0;

action.targetL0(17)=action.targetL0(17)+dL0;
action.targetL0(18)=action.targetL0(18)+dL0;
action.targetL0(19)=action.targetL0(19)+dL0;



Uhis=action.Solve();

plots.fileName='ActiveTruss.gif';
plots.activeTrussNum=[14,15,16,17,18,19];
plots.Plot_DeformedHis(Uhis);
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));



