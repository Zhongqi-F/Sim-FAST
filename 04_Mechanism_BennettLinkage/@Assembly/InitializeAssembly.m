%% This function initialize the assembly system

function InitializeAssembly(obj)

    obj.node.currentU_Mat = zeros(size(obj.node.coordinates_Mat));
    obj.node.currentExtForce_Mat = zeros(size(obj.node.coordinates_Mat));

    obj.spr.currentTheta_Vec = ...
    obj.spr.Spr_Theta(obj.node,obj.node.currentU_Mat);

    obj.spr.theta_StressFree_Vec = obj.spr.currentTheta_Vec;
end