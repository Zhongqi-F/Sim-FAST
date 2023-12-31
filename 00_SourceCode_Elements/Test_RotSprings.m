%test
clear all;
clc;
close all;

% Define the nodes
node=Elements_Nodes;
node.coordinates_Mat=[0 0 0; 1 0 0;0 1 0; 1 1 0;];


% Define the rotational springs
spr=Elements_RotSprings;
spr.rotSprIJKL_Mat=[1,2,3,4];
spr.rotSprK_Vec=[1];
spr.theta_StressFree_Vec=[pi];
spr.theta_Current_Vec=[pi];

% spr.InitializeSpr(node);


U=[0,0,0;
   0,0,0;
   0,0,0;
   0,0,0.0001;];


[T,K]=spr.SolveFK(node,U);
T
K

spr.theta_Current_Vec=[pi+0.5*pi];

[T,K]=spr.SolveFK(node,U);
T
K