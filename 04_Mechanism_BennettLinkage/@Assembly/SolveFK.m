%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=SolveFK(obj,U)
    [Tbar,Kbar]=obj.bar.SolveFK(obj.node,U);
    [Tspr,Kspr]=obj.spr.SolveFK(obj.node,U);
    
    T=Tbar+Tspr;
    K=Kbar+Kspr;
end