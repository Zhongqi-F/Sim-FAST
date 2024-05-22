%% Solve the global force and stiffness of CST

function [Tcst,Kcst]=SolveFK(obj,node,U)
    
        [Tcst]=GlobalForce(obj,node,U);
        [Kcst]=GlobalStiff(obj,node,U);
end