
classdef Assembly_Knit < handle

    properties
        % Nodes
        node

        % Bars
        bar

        % Rotational Springs
        rotSpr

        % Rotational Springs Flat
        rotSprFlat

    end

    methods
        % For a given input deformation find the global force vector and
        % the stiffness matrix
        [T,K]=SolveFK(obj,node,U)

        % Initialize the assembly
        % This will set currentU to be zero matrix
        % This will set theta_StressFree_Vec to be current theta value
        % This will set current external force to be zero vector
        InitializeAssembly(obj)

    end
end