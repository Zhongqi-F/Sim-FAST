%% Test the behavior of the bar element
clear all;
clc;
close all;

% Define the nodes
node=Elements_Nodes;
node.coordinates_Mat=[0 0 0; 1 0 0];
% Here we defined the nodal coordinates of the bar

% Define the bars
bar=Elements_Bars;

bar.barConnect_Mat=[1,2;];
bar.A_Vec=[1];
bar.E_Vec=[1];
% Here we assume that the bar is connecting node 1 and node 1
% The bar has an area of 1 and a Young's modulus of 1

bar.InitializeLengthVec(node);
% This function automaticall calculate L0 for bar element
% We could also set L0 manually

U=[0,0,0;
   0.00001,0,0;];
% Here we assume that the bar has a small extension in node 2 (x direction)

[T,K]=bar.SolveFK(node,U)
% We can see that the produced force vector and the stiffness matrix are
% all correct solution for a linear bar. 
