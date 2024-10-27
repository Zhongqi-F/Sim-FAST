%% CST elements
% This element is derived using central difference
% This is a constant strain linear elastic formulation
% The CST is constructed using 3 nodes

classdef CD_Elements_CST < handle

    properties
        % Connection information of rotational spring elements (Ns*3)
        node_ijk_mat

        % Thickness of each element (Ns*1)
        t_vec 

        % Young's modulus of each element
        E_vec;

        % Poisson's Ratio of each element
        v_vec;

        % Current Strain Energy
        energy_current_vec;

        % step for central difference
        delta=10^-8;
        
    end

    methods
        % Initialize CST elements
        Initialize(obj,node)

        % Calculate the potential of the CST elements
        PE=Potential(obj,strainMat,E,v,A,t)

        % Solve thes strain of the element
        strainMat = Solve_Strain(obj,x,Xoriginal)

        % Calculate the local force vector of the CST elements
        [Flocal]=Solve_Local_Force(obj,x,X,E,v,t)

        % Calculate the global internal force of all CST elements
        [Trsf]=Solve_Global_Force(obj,node,U)

        % Calculate the local stiffness of CST elements
        [Klocal]=Solve_Local_Stiff(obj,x,X,E,v,t)

        % Calculate the gloabl stiffness matrix of the CST elements
        [Krsf]=Solve_Global_Stiff(obj,node,U)       

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Tcst,Kcst]=Solve_FK(obj,node,U)

    end
end
