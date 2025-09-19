function [Tbar,Kbar]=Solve_FK(obj,node,U)
    
        [Tbar]=Solve_Global_Force(obj,node,U);
        [Kbar]=Solve_Global_Stiff(obj,node,U);

end