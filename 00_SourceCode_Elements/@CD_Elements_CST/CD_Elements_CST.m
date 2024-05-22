%% This function defines the CST elements. 

classdef CD_Elements_CST < handle

    properties
        % Connection information of rotational spring elements (Ns*3)
        CST_IJK_Mat

        % Thickness of each element (Ns*1)
        thickness_Vec 

        % Young's modulus of each element
        E_Vec;

        % Poisson's Ratio of each element
        v_Vec;

        % Current Strain Energy
        strainEnergy_Vec;

        % step for central difference
        delta=10^-8;
        
    end

    methods
        % Solve thes strain of the element
        strainMat = CalcStrain(obj,x,Xoriginal)

        % Calculate the potential of the CST elements
        PE=Potential(obj,strainMat,E,v,A,t)

        % Calculate the local force vector of the CST elements
        [Flocal]=LocalForce(obj,x,X,E,v,t)

        % Calculate the global internal force of all CST elements
        [Trsf]=GlobalForce(obj,node,U)

        % Calculate the local stiffness of CST elements
        [Klocal]=LocalStiff(obj,x,X,E,v,t)

        % Calculate the gloabl stiffness matrix of the CST elements
        [Krsf]=GlobalStiff(obj,node,U)

        % Initialize CST elements
        InitializeCST(obj,node)
       

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Tcst,Kcst]=SolveFK(obj,node,U)

    end
end
