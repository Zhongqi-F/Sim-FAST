function [cst_strain_mat] = ...
    Solve_CST_Strain(obj,bar_strain_Mat,trans_mat)

    % This function computes the cst strain from the bar strain
    % We directly use the factors from the trans_mat

    epsilon1_Vec = bar_strain_Mat(:,1);
    epsilon2_Vec = bar_strain_Mat(:,2);
    epsilon3_Vec = bar_strain_Mat(:,3);

    B1_vec=squeeze(trans_mat(:,2,1));
    B2_vec=squeeze(trans_mat(:,2,2));
    B3_vec=squeeze(trans_mat(:,2,3));

    C1_vec=squeeze(trans_mat(:,3,1));
    C2_vec=squeeze(trans_mat(:,3,2));
    C3_vec=squeeze(trans_mat(:,3,3));

    % We use the precomputed matrix to convert bar strain to 
    % the strain of the cst element
    epsilon_p_vec=B1_vec.*epsilon1_Vec+...
        B2_vec.*epsilon2_Vec+B3_vec.*epsilon3_Vec;
    gamma_vec=C1_vec.*epsilon1_Vec+...
        C2_vec.*epsilon2_Vec+C3_vec.*epsilon3_Vec;

    cst_strain_mat=[epsilon1_Vec,epsilon_p_vec,gamma_vec];

end