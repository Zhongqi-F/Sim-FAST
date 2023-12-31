%% This function defines the linear elastic bar elements. 
% This bar element is calculated using the central difference method, where
% we do not solve the internal force and stiffness matrix using the exact
% solution. In this element, the internal force is solved using the
% gradient of the bar potential while the stiffness is solved using the
% Hessian of the bar potential. 

classdef Elements_Bars < handle

    properties
        % Area of the bar, stored as a vector (Nb*1)
        A_Vec

        % Young's Modulus of the bar, stored as a vector (Nb*1)
        E_Vec

        % Stress Free Length of the bar, stored as a vector (Nb*1)
        L0_Vec

        % Current Strain Energy of the bar, stored as a vector (Nb*1)
        currentStrainEnergy_Vec

        % Connection information of the bar, stored as a matrix (Nb*2)
        barConnect_Mat

        % step size for evaluating gradient and hessian
        delta=10^-8;

    end

    methods
        % Potential energy of bars
        PE=Potential(obj,X1,X2,L0,E,A)
        
        % Solve the local force vector of an element
        % This is calculated as the gradient of potential
        Flocal=LocalForce(obj,X1,X2,L0,E,A)

        % Solve the local stiffness matrix of an element
        % This is calculated as the Hessian of the potential
        Klocal=LocalStiff(obj,X1,X2,L0,E,A)

        % Calculate the global force vector
        [Tbar]=GlobalForce(obj,node,U)

        % Calculate the global stiffness matrix
        [Kbar]=GlobalStiff(obj,node,U)
        
        % Initialize the original length of bar
        InitializeLengthVec(obj,node)

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tbar,Kbar]=SolveFK(obj,node,U)

        % Back calculate internal strain energy;
        CalcStrainEnergy(node,U)


    end
end
