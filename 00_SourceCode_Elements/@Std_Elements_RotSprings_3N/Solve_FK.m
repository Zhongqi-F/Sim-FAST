%% Solve the global force and stiffness of springs
% This function calcualtes the global force and stiffness of the springs

function [Trs,Krs]=Solve_FK(obj,node,U)
    
        [Trs]=Solve_Global_Force(obj,node,U);
        [Krs]=Solve_Global_Stiff(obj,node,U);
end