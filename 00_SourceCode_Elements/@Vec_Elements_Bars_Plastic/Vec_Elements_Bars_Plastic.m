%% Bar elements
% This bar element is calculated based on analytical equations
% This code is vectorized for speed
% This formulation gives a plastic response (isotropic hardening)
% The bar geometry is defined with two nodes

classdef Vec_Elements_Bars_Plastic < handle

    properties
        % Connection information of the bar, stored as a matrix (Nb*2)
        node_ij_mat

        % Area of the bar, stored as a vector (Nb*1)
        A_vec

        % Young's Modulus of the bar, stored as a vector (Nb*1)
        E_vec

        % Length of the bar, stored as a vector (Nb*1)
        L0_vec

        % Hardened Stiffness (Nb*1)
        H_vec

        % Yielding stress (Nb*1)
        sigma_y_vec

        % Current Engineering Strain of the bar, stored as a vector (Nb*1)
        strain_current_vec

        % Current Strain Energy of the bar, stored as a vector (Nb*1)
        energy_current_vec

        % Plastic Strain (Nb*1)
        strain_plastic_current_vec

        % Equivalent plastic strain (Nb*1)
        alpha_current_vec

    end

    methods
        % Calculate the strain of bars
        [Ex]=Solve_Strain(obj,node,U);

        % Calculate the local stress and stiffness of bars
        % This function defines the constitutive model of material
        [Sx,Cx]=Solve_Stress(obj,Ex)

        % Calculate the global force vector
        [Tbar]=Solve_Global_Force(obj,node,U,Sx)

        % Calculate the global stiffness matrix
        [Kbar]=Solve_Global_Stiff(obj,node,U,Sx,C)
        
        % Initialize the original length of bar
        Initialize(obj,node)

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tbar,Kbar]=Solve_FK(obj,node,U)

    end
end
