%% This function defines the linear elastic spring elements. 

classdef Elements_RotSprings < handle

    properties
        % Connection information of rotational spring elements (Ns*4)
        rotSprIJKL_Mat

        % Rotational stiffness of each element (Ns*1)
        rotSprK_Vec 

        % Stress-free angle of the spring
        theta_StressFree_Vec;

        % Current Theta
        theta_Current_Vec;

        % Current Strain Energy
        strainEnergy_Vec;

        % step for central difference
        delta=10^-8;
        
    end

    methods
        % Calculate the rotational angle of rotational springs
        theta=Theta(obj,X);

        % Calculate the potential of the rotational spring
        PE=Potential(obj,theta,theta0,K);

        % Calculate the local force vector of the rotational spring
        [Flocal]=LocalForce(obj,X,theta0,K)

        % Calculate the global internal force of all spring elements
        [Trs]=GlobalForce(obj,node,U)

        % Calculate the local stiffness of the rotational spring
        [Klocal]=LocalStiff(obj,X,theta0,K)

        % Calculate the gloabl stiffness matrix of the spring element
        [Krs]=GlobalStiff(obj,node,U)
        

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Trs,Krs]=SolveFK(obj,node,U)

        % This function initialize the springs
        InitializeSpr(obj,node)

        % Solve strain energy and current fold angle
        CalcStrainEnergy(obj,node,U)

    end
end
