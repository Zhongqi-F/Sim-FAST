%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=Solve_FK(obj,U)

    [Tbar,Kbar]=obj.bar.Solve_FK(obj.node,U);
    [Tspr,Kspr]=obj.rot_spr_3N.Solve_FK(obj.node,U);
    
    T=Tbar+Tspr;
    K=Kbar+Kspr;
end