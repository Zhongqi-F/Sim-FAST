%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Define the nodal coordinate before meshing
a=20*10^(-3);
L=50*10^(-3);
lambda=30/180*pi;
omega=60/180*pi;

% Stiffness parameters of the structure
sprStiff=0.000001;
stiffFactor=1000;
barA=10*10^(-3)*10*10^(-3); 
barE=2*10^9; % Young's modulus

% perturbation Level 
perturbationLevel=0*10^(-3);

data=zeros(30,10);

bar=Elements_Bars;
rotSpr=Elements_RotSprings;
node=Elements_Nodes;

bar.A_Vec=zeros(5,1);
bar.E_Vec=zeros(5,1);
bar.L0_Vec=zeros(5,1);


%% Bennett Linkage Geometry
% First Links
node.coordinates_Mat(1,:)=[L*cos(omega),0,0];
node.coordinates_Mat(2,:)=[L*cos(omega)+a*cos(lambda)/sin(omega),0,a*sin(lambda)];
node.coordinates_Mat(3,:)=[L*cos(omega)+a*cos(lambda)/sin(omega)+a*sin(lambda)/sin(omega),0,a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(4,:)=[L*cos(omega)+a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_Vec(1:5)=barE;
bar.A_Vec(1:5)=barA;
bar.barConnect_Mat(1,:)=[1,2];
bar.barConnect_Mat(2,:)=[2,3];
bar.barConnect_Mat(3,:)=[3,4];
bar.barConnect_Mat(4,:)=[1,4];
bar.barConnect_Mat(5,:)=[1,3];

node.coordinates_Mat(5,:)=[0,L*sin(omega),0];
node.coordinates_Mat(6,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega),a*sin(lambda)];
node.coordinates_Mat(7,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega)+a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(8,:)=[0,L*sin(omega)+a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_Vec(6:10)=barE;
bar.A_Vec(6:10)=barA;
bar.barConnect_Mat(6,:)=[5,6];
bar.barConnect_Mat(7,:)=[6,7];
bar.barConnect_Mat(8,:)=[7,8];
bar.barConnect_Mat(9,:)=[5,8];
bar.barConnect_Mat(10,:)=[5,7];

bar.E_Vec(11:14)=barE;
bar.A_Vec(11:14)=barA;
bar.barConnect_Mat(11,:)=[1,5];
bar.barConnect_Mat(12,:)=[2,6];
bar.barConnect_Mat(13,:)=[3,7];
bar.barConnect_Mat(14,:)=[4,8];

bar.E_Vec(15:18)=barE;
bar.A_Vec(15:18)=barA;
bar.barConnect_Mat(15,:)=[1,6];
bar.barConnect_Mat(16,:)=[2,7];
bar.barConnect_Mat(17,:)=[3,8];
bar.barConnect_Mat(18,:)=[4,5];

% Second Link
node.coordinates_Mat(9,:)=[0,L*sin(omega)+a*cos(lambda)/cos(omega)+a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(10,:)=[0,L*sin(omega)+a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_Vec(19:22)=barE;
bar.A_Vec(19:22)=barA;
bar.barConnect_Mat(19,:)=[6,9];
bar.barConnect_Mat(20,:)=[9,10];
bar.barConnect_Mat(21,:)=[5,10];
bar.barConnect_Mat(22,:)=[5,9];

node.coordinates_Mat(11,:)=[-L*cos(omega),0,0];
node.coordinates_Mat(12,:)=[-L*cos(omega)-a*cos(lambda)/sin(omega),0,a*sin(lambda)];
node.coordinates_Mat(13,:)=[-L*cos(omega)-a*cos(lambda)/sin(omega)-a*sin(lambda)/sin(omega),0,a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(14,:)=[-L*cos(omega)-a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_Vec(23:27)=barE;
bar.A_Vec(23:27)=barA;
bar.barConnect_Mat(23,:)=[11,12];
bar.barConnect_Mat(24,:)=[12,13];
bar.barConnect_Mat(25,:)=[13,14];
bar.barConnect_Mat(26,:)=[11,14];
bar.barConnect_Mat(27,:)=[11,13];

bar.E_Vec(28:31)=barE;
bar.A_Vec(28:31)=barA;
bar.barConnect_Mat(28,:)=[5,11];
bar.barConnect_Mat(29,:)=[6,12];
bar.barConnect_Mat(30,:)=[9,13];
bar.barConnect_Mat(31,:)=[10,14];

bar.E_Vec(32:35)=barE;
bar.A_Vec(32:35)=barA;
bar.barConnect_Mat(32,:)=[5,12];
bar.barConnect_Mat(33,:)=[6,13];
bar.barConnect_Mat(34,:)=[9,14];
bar.barConnect_Mat(35,:)=[10,11];

% Third Link
node.coordinates_Mat(15,:)=[-L*cos(omega),0,0];
node.coordinates_Mat(16,:)=[-L*cos(omega)-a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_Vec(36:39)=barE;
bar.A_Vec(36:39)=barA;
bar.barConnect_Mat(36,:)=[12,15];
bar.barConnect_Mat(37,:)=[13,16];
bar.barConnect_Mat(38,:)=[15,16];
bar.barConnect_Mat(39,:)=[12,16];

node.coordinates_Mat(17,:)=[0,-L*sin(omega),0];
node.coordinates_Mat(18,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega),a*sin(lambda)];
node.coordinates_Mat(19,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega)-a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(20,:)=[0,-L*sin(omega)-a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_Vec(40:43)=barE;
bar.A_Vec(40:43)=barA;
bar.barConnect_Mat(40,:)=[15,17];
bar.barConnect_Mat(41,:)=[12,18];
bar.barConnect_Mat(42,:)=[13,19];
bar.barConnect_Mat(43,:)=[16,20];

bar.E_Vec(44:47)=barE;
bar.A_Vec(44:47)=barA;
bar.barConnect_Mat(44,:)=[15,18];
bar.barConnect_Mat(45,:)=[12,19];
bar.barConnect_Mat(46,:)=[13,20];
bar.barConnect_Mat(47,:)=[16,17];

bar.E_Vec(48:52)=barE;
bar.A_Vec(48:52)=barA;
bar.barConnect_Mat(48,:)=[17,18];
bar.barConnect_Mat(49,:)=[18,19];
bar.barConnect_Mat(50,:)=[19,20];
bar.barConnect_Mat(51,:)=[17,20];
bar.barConnect_Mat(52,:)=[17,19];

% Forth Link
node.coordinates_Mat(21,:)=[0,-L*sin(omega)-a*cos(lambda)/cos(omega)-a*sin(lambda)/cos(omega),a*sin(lambda)-a*cos(lambda)];
node.coordinates_Mat(22,:)=[0,-L*sin(omega)-a*sin(lambda)/cos(omega),-a*cos(lambda)];

bar.E_Vec(53:56)=barE;
bar.A_Vec(53:56)=barA;
bar.barConnect_Mat(53,:)=[18,21];
bar.barConnect_Mat(54,:)=[21,22];
bar.barConnect_Mat(55,:)=[17,22];
bar.barConnect_Mat(56,:)=[17,21];

node.coordinates_Mat(23,:)=[L*cos(omega),0,0];
node.coordinates_Mat(24,:)=[L*cos(omega)+a*sin(lambda)/sin(omega),0,-a*cos(lambda)];

bar.E_Vec(57:60)=barE;
bar.A_Vec(57:60)=barA;
bar.barConnect_Mat(57,:)=[2,23];
bar.barConnect_Mat(58,:)=[3,24];
bar.barConnect_Mat(59,:)=[23,24];
bar.barConnect_Mat(60,:)=[3,23];

bar.E_Vec(61:64)=barE;
bar.A_Vec(61:64)=barA;
bar.barConnect_Mat(61,:)=[17,23];
bar.barConnect_Mat(62,:)=[2,18];
bar.barConnect_Mat(63,:)=[3,21];
bar.barConnect_Mat(64,:)=[22,24];

bar.E_Vec(65:68)=barE;
bar.A_Vec(65:68)=barA;
bar.barConnect_Mat(65,:)=[2,17];
bar.barConnect_Mat(66,:)=[2,21];
bar.barConnect_Mat(67,:)=[3,22];
bar.barConnect_Mat(68,:)=[22,23];

bar.InitializeLengthVec(node);

% Add Rotational Hinge
rotSpr.rotSprIJKL_Mat(1,:)=[2,5,6,12];
rotSpr.rotSprIJKL_Mat(2,:)=[6,12,13,18];
rotSpr.rotSprIJKL_Mat(3,:)=[12,17,18,2];
rotSpr.rotSprIJKL_Mat(4,:)=[18,2,3,6];

rotSpr.rotSprK_Vec(1,1)=stiffFactor*sprStiff;
rotSpr.rotSprK_Vec(2,1)=sprStiff;
rotSpr.rotSprK_Vec(3,1)=sprStiff;
rotSpr.rotSprK_Vec(4,1)=sprStiff;

%% Initialize assembly
assembly=Assembly_BennettLinkage();
assembly.node=node;
assembly.bar=bar;
assembly.rotSpr=rotSpr;

assembly.InitializeAssembly()

%% Plot for investigation
plot=Plot();
plot.displayRange=0.1;
plot.displayRangeRatio=1;
plot.assembly=assembly;

plot.Plot_Shape_NodeNumber()
plot.Plot_Shape_BarNumber()
plot.Plot_Shape_SprNumber()



%% Find alpha angle
% This is for comparing with the analytical solution
z1=node.coordinates_Mat(6,:)-node.coordinates_Mat(5,:);
x1=node.coordinates_Mat(6,:)-node.coordinates_Mat(2,:);
z2=node.coordinates_Mat(13,:)-node.coordinates_Mat(12,:);
x2=node.coordinates_Mat(12,:)-node.coordinates_Mat(6,:);

z1=z1/norm(z1);
x1=x1/norm(x1);
y1=cross(x1,z1);

z2=z2/norm(z2);
x2=x2/norm(x2);
y2=cross(x2,z2);

acos(dot(z1,z2))


%% Setup the loading controller
sf=Solver_NR_Folding;
sf.assembly=assembly;
sf.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;
         4,1,1,1];

sf.targetRot=rotSpr.theta_Current_Vec;
sf.targetRot(1)=sf.targetRot(1)+0.7*pi;


sf.increStep=200;
sf.tol=10^-8;
sf.iterMax=50;

Uhis=sf.Solve();

plot.displayRange=0.15;
plot.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plot.Plot_DeformedHis(Uhis)

