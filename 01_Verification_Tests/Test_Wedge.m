%% Test the behavior of the bar element
clear all;
clc;
close all;

% Define the nodes
node=Elements_Nodes;
N=10;
L=1;
H=0.1;
T=0.1;

for i=1:N
    node.coordinates_Mat=[node.coordinates_Mat;
        (i-1)*L/N, 0, 0;
        (i-1)*L/N, 0, T;
        (i-1)*L/N, H, 0;
        (i-1)*L/N, H, T;];
end


% Here we defined the nodal coordinates of the bar

% Define the bars
wedge=Elements_WedgeSolid;

wedge.barConnect_Mat=[1,2;];
wedge.A_Vec=[1];
wedge.E_Vec=[1];
% Here we assume that the bar is connecting node 1 and node 1
% The bar has an area of 1 and a Young's modulus of 1

wedge.InitializeLengthVec(node);
% This function automaticall calculate L0 for bar element
% We could also set L0 manually

U=[0,0,0;
   0.00001,0,0;];
% Here we assume that the bar has a small extension in node 2 (x direction)

[T,K]=wedge.SolveFK(node,U)
% We can see that the produced force vector and the stiffness matrix are
% all correct solution for a linear bar. 
