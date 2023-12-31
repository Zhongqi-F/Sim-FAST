%% Solve the global force and stiffness of bars
% This function calcualtes the global force and stiffness of the bars


function [Tbar,Kbar]=SolveFK(obj,node,U)
    
        [Ex]=Bar_Strain(obj,node,U);
        [Sx,C]=Bar_Cons(obj,Ex);
        [Tbar]=Bar_GlobalForce(obj,node,U,Sx);
        [Kbar]=Bar_GlobalStiffAssemble(obj,node,U,Sx,C);

        obj.currentStrain_Vec=Sx;
        obj.currentStrainEnergy_Vec=1/2*obj.E_Vec.*obj.A_Vec.*obj.L_Vec.*Ex.*Ex;

end