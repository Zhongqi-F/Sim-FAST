%% This function defines the linear elastic spring elements. 

classdef Elements_Springs < handle

    properties
        % Connection information of rotational spring elements (Ns*4)
        sprIJKL_Mat

        % Rotational stiffness of each element (Ns*1)
        sprRotK_Vec 

        % Threshold for local penertration prevention
        % The theta1 and theta2 are terms for initiating the penertration 
        % prevention. When folding angle is smaller than theta1 or bigger 
        % than theta 2, an additional forces is added. Please check the 
        % Liu and Paulino RSPA paper for details.
        theta1=0.1*pi; 
        theta2=2*pi-0.1*pi;

        % Stress-free angle of the spring
        theta_StressFree_Vec;

        % Current Theta
        currentTheta_Vec;

        % Current Strain Energy
        currentStrainEnergy_Vec;

    end

    methods
        % Calculate the rotational angle of rotational springs
        [theta]=Spr_Theta(obj,node,U);

        % Calculate the local moment and stiffness of the rotational spring
        % elements. This function gives the constitutive relationship of 
        % rotational spring elements
        [Mspr,Cspr]=Spr_Cons(obj,theta)

        % Calculate the global internal force of the spring element
        [Tspr]=Spr_GlobalForce(obj,node,U,Mspr)

        % Calculate the gloabl stiffness matrix of the spring element
        [Kspr]=Spr_GlobalStiffAssemble(obj,node,U,Mspr,Cspr)
        

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the bar elements (making use of the above four 
        % functions). 
        [Tspr,Kspr]=SolveFK(obj,node,U)

        % initialize the spring element
        InitializeSpr(obj,node)

    end
end
