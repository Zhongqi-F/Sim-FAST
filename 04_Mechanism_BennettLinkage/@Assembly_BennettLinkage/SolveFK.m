%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=SolveFK(obj,U)

    [Tbar,Kbar]=obj.bar.SolveFK(obj.node,U);
    [Trs,Krs]=obj.rotSpr.SolveFK(obj.node,U);
    
    T=Tbar+Trs;
    K=Kbar+Krs;
end