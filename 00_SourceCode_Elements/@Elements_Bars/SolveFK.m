%% Solve the global force and stiffness of bars
% This function calcualtes the global force and stiffness of the bars

function [Tbar,Kbar]=SolveFK(obj,node,U)
    
        [Tbar]=GlobalForce(obj,node,U);
        [Kbar]=GlobalStiff(obj,node,U);

end