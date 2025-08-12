
classdef Assembly_ThickOrigami < handle

    properties
        % Nodes
        node

        % CST Elements
        cst

        % Rotational Springs
        rot_spr_4N

        % Zero length springs
        zlspr

    end

    methods
        % For a given input deformation find the global force vector and
        % the stiffness matrix
        [T,K]=Solve_FK(obj,U)

        % Initialize the assembly
        % This will set currentU to be zero matrix
        % This will set theta_StressFree_Vec to be current theta value
        % This will set current external force to be zero vector
        Initialize_Assembly(obj)

        % Add the thick triangle panel
        Add_Triangle_Panel(obj,n1,n2,n3,n4,n5,n6,E,faceThickness,v)

    end
end