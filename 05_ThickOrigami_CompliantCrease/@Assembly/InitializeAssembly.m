%% This function initialize the assembly system

function InitializeAssembly(obj)

    obj.node.currentU_Mat = zeros(size(obj.node.coordinates_Mat));
    obj.node.currentExtForce_Mat = zeros(size(obj.node.coordinates_Mat));

    obj.bar.InitializeLengthVec(obj.node);
    obj.spr.InitializeSpr(obj.node);

    obj.wedge.InitializeL0ref(obj.node);
end