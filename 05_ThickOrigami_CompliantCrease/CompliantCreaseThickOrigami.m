%% Initialize the solver
clear all;clc;close all;


%% Define the geometry of thick origami
% This section of code is used to generate geometry of thick origami

% Define the nodal coordinate before meshing
L=50*10^(-3);
alpha=30/180*pi;
t=5*10^(-3);
gap=6*10^-3;


% Stiffness parameters of the structure
sprStiff=10^-6;
stiffFactor=10000;
barA=5*10^(-3)*30*10^(-3);
barE=2*10^9; % Young's modulus
wedgeK=10000; % stiffness of the wedge

% perturbation Level 
perturbationLevel=0.05*10^(-3);

data=zeros(30,10);

bar=Elements_Bars;
spr=Elements_Springs;
node=Elements_Nodes;
wedge=Elements_WedgeSolid;

%% thick origami geometry

% First panel
node.coordinates_Mat(1,:)=[gap/2/tan(22.5/180*pi),gap/2,0];
node.coordinates_Mat(2,:)=[L,gap/2,0];
node.coordinates_Mat(3,:)=[L,L-gap/sqrt(2),0];

node.coordinates_Mat(4,:)=[gap/2/tan(22.5/180*pi),gap/2,t];
node.coordinates_Mat(5,:)=[L,gap/2,t];
node.coordinates_Mat(6,:)=[L,L-gap/sqrt(2),t];

% Second panel
node.coordinates_Mat(7,:)=[gap/2/tan(22.5/180*pi)-gap/sqrt(2),gap/2+gap/sqrt(2),0];
node.coordinates_Mat(8,:)=[L-gap/sqrt(2),L,0];
node.coordinates_Mat(9,:)=[-L+gap/sqrt(2),L,0];
node.coordinates_Mat(10,:)=[-gap/2/tan(22.5/180*pi)+gap/sqrt(2),gap/2+gap/sqrt(2),0];

node.coordinates_Mat(11,:)=[gap/2/tan(22.5/180*pi)-gap/sqrt(2),gap/2+gap/sqrt(2),t];
node.coordinates_Mat(12,:)=[L-gap/sqrt(2),L,t];
node.coordinates_Mat(13,:)=[-L+gap/sqrt(2),L,t];
node.coordinates_Mat(14,:)=[-gap/2/tan(22.5/180*pi)+gap/sqrt(2),gap/2+gap/sqrt(2),t];

% Third panel
node.coordinates_Mat(15,:)=[-gap/2/tan(22.5/180*pi),gap/2,0];
node.coordinates_Mat(16,:)=[-L,gap/2,0];
node.coordinates_Mat(17,:)=[-L,L-gap/sqrt(2),0];

node.coordinates_Mat(18,:)=[-gap/2/tan(22.5/180*pi),gap/2,t];
node.coordinates_Mat(19,:)=[-L,gap/2,t];
node.coordinates_Mat(20,:)=[-L,L-gap/sqrt(2),t];

% Forth panel
node.coordinates_Mat(21,:)=[-gap/2/tan(22.5/180*pi),-gap/2,0];
node.coordinates_Mat(22,:)=[-L,-gap/2,0];
node.coordinates_Mat(23,:)=[-L,-L+gap/sqrt(2),0];

node.coordinates_Mat(24,:)=[-gap/2/tan(22.5/180*pi),-gap/2,t];
node.coordinates_Mat(25,:)=[-L,-gap/2,t];
node.coordinates_Mat(26,:)=[-L,-L+gap/sqrt(2),t];

% Fifth panel
node.coordinates_Mat(27,:)=[gap/2/tan(22.5/180*pi)-gap/sqrt(2),-gap/2-gap/sqrt(2),0];
node.coordinates_Mat(28,:)=[L-gap/sqrt(2),-L,0];
node.coordinates_Mat(29,:)=[-L+gap/sqrt(2),-L,0];
node.coordinates_Mat(30,:)=[-gap/2/tan(22.5/180*pi)+gap/sqrt(2),-gap/2-gap/sqrt(2),0];

node.coordinates_Mat(31,:)=[gap/2/tan(22.5/180*pi)-gap/sqrt(2),-gap/2-gap/sqrt(2),t];
node.coordinates_Mat(32,:)=[L-gap/sqrt(2),-L,t];
node.coordinates_Mat(33,:)=[-L+gap/sqrt(2),-L,t];
node.coordinates_Mat(34,:)=[-gap/2/tan(22.5/180*pi)+gap/sqrt(2),-gap/2-gap/sqrt(2),t];

% Sixth panel
node.coordinates_Mat(35,:)=[gap/2/tan(22.5/180*pi),-gap/2,0];
node.coordinates_Mat(36,:)=[L,-gap/2,0];
node.coordinates_Mat(37,:)=[L,-L+gap/sqrt(2),0];

node.coordinates_Mat(38,:)=[gap/2/tan(22.5/180*pi),-gap/2,t];
node.coordinates_Mat(39,:)=[L,-gap/2,t];
node.coordinates_Mat(40,:)=[L,-L+gap/sqrt(2),t];

% Set up information of the wedge
wedge.wedgeConnect_Mat=[1,2,3,4,5,6;
                        7,8,9,11,12,13;
                        7,9,10,11,13,14;
                        15,16,17,18,19,20;
                        21,22,23,24,25,26;
                        27,28,29,31,32,33;
                        27,29,30,31,33,34;
                        35,36,37,38,39,40];
wedge.Kwedge_Vec=[wedgeK;
                  wedgeK;
                  wedgeK;
                  wedgeK;
                  wedgeK;
                  wedgeK;
                  wedgeK;
                  wedgeK];
wedge.InitializeL0ref(node);

% Setup bars and rotational springs for connection
bar.barConnect_Mat=zeros(1,2);
spr.sprIJKL_Mat=zeros(1,4);

bar.barConnect_Mat(1,:)=[4,11];
bar.barConnect_Mat(2,:)=[6,12];
bar.barConnect_Mat(3,:)=[4,12];

spr.sprIJKL_Mat(1,:)=[12,4,6,5];
spr.sprIJKL_Mat(2,:)=[11,4,12,6];
spr.sprIJKL_Mat(3,:)=[13,11,12,4];

bar.barConnect_Mat(4,:)=[14,18];
bar.barConnect_Mat(5,:)=[13,20];
bar.barConnect_Mat(6,:)=[14,20];

spr.sprIJKL_Mat(4,:)=[11,13,14,20];
spr.sprIJKL_Mat(5,:)=[13,14,20,18];
spr.sprIJKL_Mat(6,:)=[14,20,18,19];

bar.barConnect_Mat(7,:)=[15,21];
bar.barConnect_Mat(8,:)=[16,22];
bar.barConnect_Mat(9,:)=[15,22];

spr.sprIJKL_Mat(7,:)=[17,16,15,22];
spr.sprIJKL_Mat(8,:)=[16,22,15,21];
spr.sprIJKL_Mat(9,:)=[15,21,22,23];

bar.barConnect_Mat(10,:)=[24,34];
bar.barConnect_Mat(11,:)=[26,33];
bar.barConnect_Mat(12,:)=[24,33];

spr.sprIJKL_Mat(10,:)=[25,26,24,33];
spr.sprIJKL_Mat(11,:)=[26,33,24,34];
spr.sprIJKL_Mat(12,:)=[24,34,33,32];

bar.barConnect_Mat(13,:)=[31,38];
bar.barConnect_Mat(14,:)=[32,40];
bar.barConnect_Mat(15,:)=[31,40];

spr.sprIJKL_Mat(13,:)=[33,32,31,40];
spr.sprIJKL_Mat(14,:)=[32,40,31,38];
spr.sprIJKL_Mat(15,:)=[31,38,40,39];

bar.barConnect_Mat(16,:)=[1,35];
bar.barConnect_Mat(17,:)=[2,36];
bar.barConnect_Mat(18,:)=[1,36];

spr.sprIJKL_Mat(16,:)=[37,36,35,1];
spr.sprIJKL_Mat(17,:)=[35,1,36,2];
spr.sprIJKL_Mat(18,:)=[36,1,2,3];

bar.A_Vec=barA*ones(18,1);
bar.E_Vec=barE*ones(18,1);

spr.theta_StressFree_Vec=pi*ones(18,1);
spr.sprRotK_Vec=sprStiff*ones(18,1);
spr.sprRotK_Vec([2,5,8,11,14,17])=spr.sprRotK_Vec([2,5,8,11,14,17])*stiffFactor;



%% Initialize assembly

assembly=Assembly();
assembly.wedge=wedge;
assembly.node=node;
assembly.bar=bar;
assembly.spr=spr;
assembly.InitializeAssembly();



%% Plot for investigation
plots=Plot();
plots.displayRange=0.1;
plots.displayRangeRatio=1;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber();
plots.Plot_Shape_BarNumber();
plots.Plot_Shape_SprNumber();



%% Include imperfection at nodal coordinates

% v_perturb1=rand(3,1);
% v_perturb1=(v_perturb1-0.5);
% v_perturb1=v_perturb1/norm(v_perturb1);
% v_perturb1=v_perturb1*perturbationLevel;
% 
% v_perturb2=rand(3,1);
% v_perturb2=(v_perturb2-0.5);
% v_perturb2=v_perturb2/norm(v_perturb2);
% v_perturb2=v_perturb2*perturbationLevel;
% 
% node.coordinates_Mat(12,:)=node.coordinates_Mat(12,:)+v_perturb1';
% node.coordinates_Mat(13,:)=node.coordinates_Mat(13,:)+v_perturb2';


%% check stiffness
[F,K]=assembly.SolveFK(zeros(40,3));
[V,U]=eigs(K,10,'smallestabs')

%% Setup the loading controller
sf=Solver_NR_Folding;
sf.assembly=assembly;
sf.supp=[31,1,1,1;
         32,1,1,1;
         33,1,1,1;
         34,1,1,1;];

sf.targetRot=spr.currentTheta_Vec;
originalRot=spr.currentTheta_Vec;

targetRate=0.4;

sf.targetRot(1)=sf.targetRot(1)-targetRate*pi;
sf.targetRot(3)=sf.targetRot(1)-targetRate*pi;

sf.targetRot(4)=sf.targetRot(4)-targetRate*pi;
sf.targetRot(6)=sf.targetRot(6)-targetRate*pi;

sf.targetRot(7)=sf.targetRot(7)+targetRate*pi;
sf.targetRot(9)=sf.targetRot(9)-targetRate*pi;


sf.targetRot(10)=sf.targetRot(10)-targetRate*pi;
sf.targetRot(12)=sf.targetRot(12)+targetRate*pi;

sf.targetRot(13)=sf.targetRot(13)-targetRate*pi;
sf.targetRot(15)=sf.targetRot(15)+targetRate*pi;

sf.targetRot(16)=sf.targetRot(16)+targetRate*pi;
sf.targetRot(18)=sf.targetRot(18)-targetRate*pi;


sf.increStep=100;
sf.tol=10^-7;
sf.iterMax=100;

Uhis=sf.Solve();

plots.displayRange=0.15;
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
% plot.Plot_DeformedHis(Uhis)


%% Solve the folding angle 
angleHis=zeros(sf.increStep,18);
momentHis=zeros(sf.increStep,18);
stressFreeAngleHis=zeros(sf.increStep,1);

for i=1:sf.increStep
    angle=spr.Spr_Theta(node,squeeze(Uhis(i,:,:)));
    angleHis(i,:)=angle;

    rot1TargetFold=i/sf.increStep*sf.targetRot(1)+(1-i/sf.increStep)*originalRot(1);
    spr.theta_StressFree_Vec(1)=rot1TargetFold;
    stressFreeAngleHis(i)=rot1TargetFold;

    [moment,Krot]=spr.Spr_Cons(angle);
    momentHis(i,:)=moment;
end

Output=zeros(sf.increStep,2);
Output(:,1)=angleHis(:,1);
Output(:,2)=-momentHis(:,1);

figure
plot(stressFreeAngleHis,Output(:,2));