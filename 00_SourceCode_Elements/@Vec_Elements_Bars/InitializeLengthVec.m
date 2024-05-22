
function InitializeLengthVec(obj,node)

    obj.L0_Vec=zeros(size(obj.A_Vec));

    for i=1:length(obj.A_Vec)
        node1=node.coordinates_Mat(obj.barConnect_Mat(i,1),:);
        node2=node.coordinates_Mat(obj.barConnect_Mat(i,2),:);

        obj.L0_Vec(i)=norm(node1-node2);        
    end

end
