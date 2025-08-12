clear all
close all
clc

tic

%% Define Geometry
N=4;
L=0.1;
% L=0.2;
% L=0.05;
W=0.1;

node=Elements_Nodes;

for i=0:N
    node.coordinates_mat=[node.coordinates_mat;
        i*L   W   0;
        i*L   0   0];
end

for i=0:N
    node.coordinates_mat=[node.coordinates_mat;
        i*L   3*W   0;
        i*L   2*W   0];
end



%% Define assembly
assembly=Assembly_cstvsbar;
cst=Vec_Elements_CST;
bar=Vec_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;

%% Define Plotting Functions
plots=Plot_cstvsbar;
plots.assembly=assembly;
plots.displayRange=[-0.1; 5*L; -0.1; 0.5; -0.3; 0.3];

plots.viewAngle1=90;
plots.viewAngle2=90;

plots.Plot_Shape_NodeNumber;


%% Define Triangle

cst.node_ijk_mat=[
    1 2 3;
    2 3 4;
    3 4 6;
    3 6 5;
    5 6 7;
    6 7 8;
    7 8 10;
    7 9 10;];

% thickness of the triangle
t=0.001;
E=10^9;
cst.t_vec=t*ones((N)*2,1); % 1mm thickness
cst.E_vec=E*ones((N)*2,1); % 10^8 Young's Modulus
cst.v_vec=0.2*ones((N)*2,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;
plots.Plot_Shape_NodeNumber;

%% Set up bar element
bar.node_ij_mat=[
    11 12;
    11 13;
    12 13;
    12 14;
    13 14;
    13 16;
    14 16;
    13 15;
    15 16;
    16 17;
    15 17;
    16 18;
    17 18;
    17 20; 
    18 20;
    17 19;
    19 20;];


bar.E_vec=E*ones(17,1);
bar.A_vec=ones(17,1)*t*2*(W*4*L)/(1-0.2)/(8*L+5*W+4*sqrt(L^2+W^2));

plots.Plot_Shape_BarNumber;
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    [1:20]', zeros(20,1), zeros(20,1), ones(20,1)
];
nr.supp(1,:)=[1 1 1 1];
nr.supp(2,:)=[2 1 1 1];
nr.supp(11,:)=[11 1 1 1];
nr.supp(12,:)=[12 1 1 1];

force=100;
step=5;

Uhis=[];
for k=1:step
    nr.load=[9,force*k/step,0,0,;
             10,force*k/step,0,0,;
             19,force*k/step,0,0,;
             20,force*k/step,0,0,;];
    
    nr.increStep=1;
    nr.iterMax=10;
    nr.tol=1*10^-4;

    Uhis(k,:,:)=squeeze(nr.Solve());
end



toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='Verification_cstvsbar_extension.gif';
plots.Plot_DeformedHis(Uhis)



%% F-disp curve
UhisBar=squeeze(Uhis(:,[19,20],1));
UhisBar=mean(UhisBar,2);
UhisBar=[0;UhisBar];

UhisCst=squeeze(Uhis(:,[9,10],1));
UhisCst=mean(UhisCst,2);
UhisCst=[0;UhisCst];

forceHis=(1:step)*force*2/step;
forceHis=forceHis';
forceHis=[0;forceHis];

% theoretical solution
Krib=E*t*W/(4*L);
dth=0.001;
fth=dth*Krib;

figure
hold on
plot(UhisBar,forceHis);
plot(UhisCst,forceHis);
plot([0,dth],[0,fth]);

xlabel('displacement (m)')
ylabel('force (N)')
legend('simulation-bar','simulation-cst','theory')

