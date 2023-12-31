%% Solve the global force and stiffness of springs
% This function calcualtes the global force and stiffness of the springs


function [Tspr,Kspr]=SolveFK(obj,node,U)
    
        [theta]=Spr_Theta(obj,node,U);
        [Mspr,Cspr]=Spr_Cons(obj,theta);
        [Tspr]=Spr_GlobalForce(obj,node,U,Mspr);
        [Kspr]=Spr_GlobalStiffAssemble(obj,node,U,Mspr,Cspr);

        obj.currentTheta_Vec=theta;
        obj.currentStrainEnergy_Vec=1/2*obj.sprRotK_Vec.*...
            (theta-obj.theta_StressFree_Vec).*...
            (theta-obj.theta_StressFree_Vec);


end