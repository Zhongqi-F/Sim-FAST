%% 3-Node rotational spring elements. 
% This element is derived using central difference method
% This is a linear elastic rotational spring
% This rotational spring is defined using three nodes

classdef CD_Elements_RotSprings_3N < handle

    properties
        % Connection information of rotational spring elements (Ns*3)
        node_ijk_mat

        % Rotational stiffness of each element (Ns*1)
        rot_spr_K_vec 

        % Stress-free angle of the spring
        theta_stress_free_vec;

        % Current Theta
        theta_current_vec;

        % Current Strain Energy
        energy_current_vec;

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
        [Trsf]=GlobalForce(obj,node,U)

        % Calculate the local stiffness of the rotational spring
        [Klocal]=LocalStiff(obj,X,theta0,K)

        % Calculate the gloabl stiffness matrix of the spring element
        [Krsf]=GlobalStiff(obj,node,U)
        

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Trsf,Krsf]=SolveFK(obj,node,U)

        % This function initialize the springs
        InitializeSpr(obj,node)

        % Solve strain energy and current fold angle
        CalcStrainEnergy(obj,node,U)

    end
end
