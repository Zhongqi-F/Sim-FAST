function [Tcst,Kcst]=Solve_FK(obj,node,U)
    
        [Tcst]=Solve_Global_Force(obj,node,U);
        [Kcst]=Solve_Global_Stiff(obj,node,U);
end