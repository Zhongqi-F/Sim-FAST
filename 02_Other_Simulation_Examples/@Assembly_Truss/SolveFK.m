%% Find the global internal force and stiffness matrix of the assembly

function [T,K]=SolveFK(obj,U)

    [Tbar,Kbar]=obj.bar.Solve_FK(obj.node,U);
    
    T=Tbar;
    K=Kbar;
end