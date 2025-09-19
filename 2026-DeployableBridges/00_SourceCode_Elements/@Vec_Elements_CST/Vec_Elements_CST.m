%% CST elements
% This element is derived using analytical equations
% This code is vectorized for speed 
% This is a constant strain linear elastic formulation
% The CST is constructed using 3 nodes

classdef Vec_Elements_CST < handle

    properties
        % Connection information of rotational spring elements (Ncst*3)
        node_ijk_mat

        % Thickness of each element (Ncst*1)
        t_vec 

        % Young's modulus of each element (Ncst*1)
        E_vec;

        % Poisson's Ratio of each element (Ncst*1)
        v_vec;

        % Triangle Area of each element (Ncst*1)
        A_vec;

        % Original Length of each side (Ncst*3)
        L_mat;

        % Current Strain Energy (Ncst*1)
        energy_current_vec;
        
    end

    methods

        % Initialize CST elements
        % This include solve for the A_Vec, alpha_Vec,
        % L_Mat, and the Trans_Mat.
        Initialize(obj,node)

        % Solve the bar strain of the element
        [bar_strain_mat,l_mat,x_reshape,trans_mat] = ...
            Solve_Bar_Strain(obj,U,X0)

        % Solve thes strain of the element
        [cst_strain_Mat] = ...
            Solve_CST_Strain(obj,bar_strain_mat,trans_mat)

        % Solve the derivatives for dedx and d2edx2
        [dedx,d2edx2] = Solve_Derivatives(obj,x_reshape,l_mat)

        % Calculate the global internal force of all CST elements
        [Trsf]=Solve_Global_Force(obj,U,dedx,cst_strain_mat,trans_mat)

        % Calculate the gloabl stiffness matrix of the CST elements
        [Krsf]=Solve_Global_Stiff(obj,U,dedx,d2edx2,cst_strain_mat,trans_mat)

        % This function is the main function we use to interact with the
        % solver. We use this function to compute the global forces and 
        % stiffness of the rotational spring elements (making use of the 
        % above four functions). 
        [Tcst,Kcst]=Solve_FK(obj,node,U)

    end
end
