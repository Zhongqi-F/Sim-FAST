%% Constitutive relationship of bars
% This function calculate the stress and material stiffness once the
% strain and material properties are given

function  [Sx,Cx]= Bar_Cons(obj,Ex)

    Sx=obj.E_Vec.*Ex;
    Cx=obj.E_Vec;

end