
classdef Assembly_ThickOrigami < handle

    properties
        % Nodes
        node

        % Bars
        bar

        % Rotational Springs
        rotSpr

        % Wedges
        wedge

    end

    methods
        % For a given input deformation find the global force vector and
        % the stiffness matrix
        [T,K]=SolveFK(obj,node,U)

        % Add a compliant crease to the specified location
        AddCompliantCrease(obj,LeftNode,RightNode,LeftAnchor,RightAnchor,t,W,E,G)

        % Initialize the assembly
        % This will set currentU to be zero matrix
        % This will set theta_StressFree_Vec to be current theta value
        % This will set current external force to be zero vector
        InitializeAssembly(obj)

    end
end