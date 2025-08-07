%% Bar elements 
% This bar element is calculated using the central difference method
% This formulation gives a linear elastic response
% The bar geometry is defined with two nodes

classdef Std_Elements_Bars < handle

    properties
        % Connection information of the bar, stored as a matrix (Nb*2)
        node_ij_mat

        % Area of the bar, stored as a vector (Nb*1)
        A_vec

        % Young's Modulus of the bar, stored as a vector (Nb*1)
        E_vec

        % Stress Free Length of the bar, stored as a vector (Nb*1)
        L0_vec

        % step size for evaluating gradient and hessian
        delta=10^-8;

    end

    methods
        % Initialize the original length of bar
        Initialize(obj,node)
        
        % Solve the local force vector of an element
        % This is calculated as the gradient of potential
        Flocal=Solve_Local_Force(obj,X1,X2,L0,E,A)

        % Solve the local stiffness matrix of an element
        % This is calculated as the Hessian of the potential
        Klocal=Solve_Local_Stiff(obj,X1,X2,L0,E,A)

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
