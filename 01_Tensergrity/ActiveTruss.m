%test
clear all;
clc;
close all;

% Define the nodes
node=Elements_Nodes;
node.coordinates_Mat=[0 0 0; 1 0 0;0 1 0; 0 0 1;];

% Define the bars
bar=Elements_Bars_CentralDifference;

bar.barConnect_Mat=[1,4;
                    2,4;
                    3,4;];

bar.A_Vec=ones(3,1);
bar.E_Vec=ones(3,1);

bar.InitializeLengthVec(node);


assembly=Assembly_Truss;
assembly.bar=bar;
assembly.node=node;

assembly.InitializeAssembly();


plot=Plot();
plot.displayRangeRatio=1;
plot.assembly=assembly;

plot.Plot_Shape_NodeNumber()
plot.Plot_Shape_BarNumber()


action=Solver_NR_TrussAction;

action.assembly=assembly;
action.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;];

action.targetL0=bar.L0_Vec;
action.targetL0(3)=action.targetL0(3)+0.5;

Uhis=action.Solve();

plot.Plot_DeformedHis(Uhis);
plot.Plot_DeformedShape(squeeze(Uhis(end,:,:)));



