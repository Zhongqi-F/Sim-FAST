%% Solve the global force and stiffness of springs
% This function calcualtes the global force and stiffness of the springs


function [Trs,Krs]=SolveFK(obj,node,U)
    
        [Trs]=GlobalForce(obj,node,U);
        [Krs]=GlobalStiff(obj,node,U);
end