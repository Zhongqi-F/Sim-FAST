clear all
close all
clc

tic

%% Define Geometry
N=1;
% L=0.1;
L=0.2;
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
plots.displayRange=[-0.1; (N+1)*L; -0.1; 0.5; -0.3; 0.3];

plots.viewAngle1=90;
plots.viewAngle2=90;

plots.Plot_Shape_NodeNumber;


%% Define Triangle

cst.node_ijk_mat=[
    1 2 3;
    2 3 4;];

% thickness of the triangle
t=0.001;
E=10^9;
v=0.2;
cst.t_vec=t*ones((N)*2,1); % 1mm thickness
cst.E_vec=E*ones((N)*2,1); % 10^8 Young's Modulus
cst.v_vec=v*ones((N)*2,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;
plots.Plot_Shape_NodeNumber;

%% Set up bar element
bar.node_ij_mat=[
    5 6;
    5 7;
    6 8;
    7 8;
    5 8;
];


bar.E_vec=E*ones(5,1);
bar.A_vec=ones(5,1)*t*2*(W*L)/(1-v)/(2*L+2*W+1*sqrt(L^2+W^2));

plots.Plot_Shape_BarNumber;
assembly.Initialize_Assembly;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[
    [1:8]', zeros(8,1), zeros(8,1), ones(8,1)
];
nr.supp(1,:)=[1 1 1 1];
nr.supp(2,:)=[2 1 1 1];
nr.supp(5,:)=[5 1 1 1];
nr.supp(6,:)=[6 1 1 1];

nr.supp(3,:)=[3 1 0 1];
nr.supp(4,:)=[4 1 0 1];
nr.supp(7,:)=[7 1 0 1];
nr.supp(8,:)=[8 1 0 1];

force=10;
step=5;

Uhis=[];
for k=1:step
    nr.load=[3,0,force*k/step,0,;
             4,0,force*k/step,0,;
             7,0,force*k/step,0,;
             8,0,force*k/step,0,;];
    
    nr.increStep=1;
    nr.iterMax=10;
    nr.tol=1*10^-4;

    Uhis(k,:,:)=squeeze(nr.Solve());
end



toc

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='Verification_cstvsbar_shear.gif';
plots.Plot_DeformedHis(Uhis)



%% F-disp curve
UhisBar=squeeze(Uhis(:,[7,8],2));
UhisBar=mean(UhisBar,2);
UhisBar=[0;UhisBar];

UhisCst=squeeze(Uhis(:,[3,4],2));
UhisCst=mean(UhisCst,2);
UhisCst=[0;UhisCst];

forceHis=(1:step)*force*2/step;
forceHis=forceHis';
forceHis=[0;forceHis];

% theoretical solution
G=E/(2*(1+v));
Krib=G*t*W/(L);
dth=0.0001;
fth=dth*Krib;

figure
hold on
plot(UhisBar,forceHis);
plot(UhisCst,forceHis);
plot([0,dth],[0,fth]);

xlabel('displacement (m)')
ylabel('force (N)')
legend('simulation-bar','simulation-cst','theory')

