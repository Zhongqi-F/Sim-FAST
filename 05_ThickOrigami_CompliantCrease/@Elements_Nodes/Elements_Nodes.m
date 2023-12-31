classdef Elements_Nodes < handle

    properties
        % Nodal corridnates stored as a matrix (N*3)
        % Total number of node (N)
        coordinates_Mat

        % Mass of nodes, stored as a vector (N*1)
        mass_Vec

        % Current deformation, stored as a Matrix (N*3)
        currentU_Mat

        % Current external force (N*3)
        currentExtForce_Mat

    end

    methods
        



    end
end