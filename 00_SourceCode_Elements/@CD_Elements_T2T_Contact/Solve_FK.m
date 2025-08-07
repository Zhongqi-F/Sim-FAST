function [Tcontact,Kcontact]=Solve_FK(obj,node,U)

        [Tcontact]=obj.Solve_Global_Force(node,U);
        [Kcontact]=obj.Solve_Global_Stiff(node,U);

end