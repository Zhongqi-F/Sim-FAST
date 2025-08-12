clear all
close all
clc

tic

%% Define Geometry
N=10000;
L=0.1;
W=0.1;

node=Elements_Nodes;

for i=0:N
    node.coordinates_mat=[node.coordinates_mat;
        i*L   W   0;
        i*L   0   0];
end


%% Define assembly
assembly=Assembly_bar;
bar=Vec_Elements_Bars;
assembly.node=node;
assembly.bar=bar;

%% Define Plotting Functions
plots=Plot_bar;
plots.assembly=assembly;
plots.displayRange=[-0.1; 5*L; -0.1; 0.5; -0.3; 0.3];

plots.viewAngle1=90;
plots.viewAngle2=90;

plots.Plot_Shape_NodeNumber;


%% Set up bar element

for i=1:N
    bar.node_ij_mat=[bar.node_ij_mat;
        2*(i-1)+1 2*(i-1)+2
        2*(i-1)+2 2*(i-1)+3
        2*(i-1)+2 2*(i-1)+4
        2*(i-1)+1 2*(i-1)+3];
end

bar.node_ij_mat=[bar.node_ij_mat;
    2*N+1    2*N+2];

t=0.001;
E=10^9;
bar.E_vec=E*ones(4*N+1,1);
bar.A_vec=ones(4*N+1,1)*t*2*(W*4*L)/(1-0.2)/(8*L+5*W+4*sqrt(L^2+W^2));

plots.Plot_Shape_BarNumber;
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[[1:2*(N+1)]', zeros(2*(N+1),1), zeros(2*(N+1),1), ones(2*(N+1),1)];
nr.supp(1,:)=[1 1 1 1];
nr.supp(2,:)=[2 1 1 1];

force=100;
step=5;

Uhis=[];
for k=1:step
    nr.load=[2*(N+1),force*k/step,0,0;
             2*(N+1)-1,force*k/step,0,0];
    
    nr.increStep=1;
    nr.iterMax=10;
    nr.tol=1*10^-4;

    Uhis(k,:,:)=squeeze(nr.Solve());
end

toc




