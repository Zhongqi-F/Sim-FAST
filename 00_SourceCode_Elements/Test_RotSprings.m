%test
clear all;
clc;
close all;

% Define the nodes
node=Elements_Nodes;
node.coordinates_Mat=[0 0 0; 1 0 0;0 1 0; 0 -1 0;];

% Define the bars
bar=Elements_Bars_CentralDifference;
% bar=Elements_Bars;

bar.barConnect_Mat=[1,2;
                    1,3;
                    1,4;
                    2,3;
                    2,4;];
bar.A_Vec=ones(5,1);
bar.E_Vec=ones(5,1);

bar.InitializeLengthVec(node);


% Define the rotational springs
spr=Elements_SpringsNonVectorized;
spr.sprIJKL_Mat=[4,1,2,3];
spr.sprRotK_Vec=[1];
spr.theta_StressFree_Vec=[pi];


U=[0,0,0;
   0.0001,0,0;
   0,0,0;
   0,0,0;];

assembly=Assembly;
assembly.bar=bar;
assembly.spr=spr;
assembly.node=node;

assembly.InitializeAssembly();


[T,K]=assembly.SolveFK(U)


plot=Plot();
plot.displayRangeRatio=1;
plot.assembly=assembly;

plot.Plot_Shape_NodeNumber()
plot.Plot_Shape_BarNumber()
plot.Plot_Shape_SprNumber()



nr=Solver_NR_Loading;

nr.assembly=assembly;
nr.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;];
nr.load=[4,0,0,0.01];

Uhis=nr.Solve();

plot.Plot_DeformedHis(Uhis);
plot.Plot_DeformedShape(squeeze(Uhis(end,:,:)));


sf=Solver_NR_Folding;

sf.assembly=assembly;
sf.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;];


sf.targetRot=[pi+0.5*pi];

Uhis=sf.Solve();

% plot.Plot_DeformedHis(Uhis);
plot.Plot_DeformedShape(squeeze(Uhis(end,:,:)));




