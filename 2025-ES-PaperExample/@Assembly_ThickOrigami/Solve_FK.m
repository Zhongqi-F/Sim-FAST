%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=Solve_FK(obj,U)

    [Tcst,Kcst]=obj.cst.Solve_FK(obj.node,U);
    [Trs,Krs]=obj.rot_spr_4N.Solve_FK(obj.node,U);
    [Tzlspr,Kzlspr]=obj.zlspr.Solve_FK(obj.node,U);
    
    T=Tcst+Trs+Tzlspr;
    K=Kcst+Krs+Kzlspr;
end