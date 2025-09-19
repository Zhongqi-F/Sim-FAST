function [Tcst,Kcst]=Solve_FK(obj,node,U)

        [bar_strain_mat,l_mat,x_reshape,trans_mat]=...
            obj.Solve_Bar_Strain(U,node.coordinates_mat);

        [cst_strain_mat] = obj.Solve_CST_Strain(bar_strain_mat,trans_mat);

        [dedx,d2edx2] = obj.Solve_Derivatives(x_reshape,l_mat);
    
        [Tcst]=obj.Solve_Global_Force(U,dedx,cst_strain_mat,trans_mat);

        [Kcst]=obj.Solve_Global_Stiff(U,dedx,d2edx2,cst_strain_mat,trans_mat);
        
end