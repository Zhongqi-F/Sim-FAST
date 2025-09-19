%% Zero Length Spring Elements 
% This element is calculated based on analytical equations
% This code is vectorized for speed
% This formulation gives a linear elastic response
% The bar geometry is defined with two nodes

classdef Vec_Elements_Zero_L_Spring < handle

    properties
        % Connection information of the zero length spring, stored as a matrix (Nb*2)
        node_ij_mat

        % Stiffness of the spring, stored as a vector (Nb*1)
        k_vec

        % Current Engineering Strain of the bar, stored as a vector (Nb*1)
        extension_current_vec

        % Current Strain Energy of the bar, stored as a vector (Nb*1)
        energy_current_vec

    end

    methods
        % Initialize the original length of bar
        Initialize(obj,node)

        % Calculate the global force vector
        [Tbar]=Solve_Global_Force(obj,node,U)

        % Calculate the global stiffness matrix
        [Kbar]=Solve_Global_Stiff(obj,node,U)
        
        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tbar,Kbar]=Solve_FK(obj,node,U)

    end
end
