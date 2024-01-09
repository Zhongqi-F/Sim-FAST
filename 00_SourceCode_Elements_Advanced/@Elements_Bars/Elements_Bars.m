%% This function defines the linear elastic bar elements. 

classdef Elements_Bars < handle

    properties
        % Area of the bar, stored as a vector (Nb*1)
        A_Vec

        % Young's Modulus of the bar, stored as a vector (Nb*1)
        E_Vec

        % Length of the bar, stored as a vector (Nb*1)
        L0_Vec

        % Current Engineering Strain of the bar, stored as a vector (Nb*1)
        currentStrain_Vec

        % Current Strain Energy of the bar, stored as a vector (Nb*1)
        currentStrainEnergy_Vec

        % Connection information of the bar, stored as a matrix (Nb*2)
        barConnect_Mat

    end

    methods
        % Calculate the strain of bars
        [Ex]=Bar_Strain(obj,node,U);

        % Calculate the local stress and stiffness of bars
        % This function defines the constitutive model of material
        [Sx,Cx]=Bar_Cons(obj,Ex)

        % Calculate the global force vector
        [Tbar]=Bar_GlobalForce(obj,node,U,Sx)

        % Calculate the global stiffness matrix
        [Kbar]=Bar_GlobalStiffAssemble(obj,node,U,Sx,C)
        
        % Initialize the original length of bar
        InitializeLengthVec(obj,node)

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tbar,Kbar]=SolveFK(obj,node,U)

    end
end
