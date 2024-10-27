function [Tzlspr,Kzlspr]=Solve_FK(obj,node,U)
    
        [Tzlspr]=Solve_Global_Force(obj,node,U);
        [Kzlspr]=Solve_Global_Stiff(obj,node,U);

end