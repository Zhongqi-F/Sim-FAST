%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=SolveFK(obj,U)

    [Tcst,Kcst]=obj.cst.SolveFK(obj.node,U);
    [Trs,Krs]=obj.rotSpr.SolveFK(obj.node,U);
    
    T=Tcst+Trs;
    K=Kcst+Krs;
end