
function InitializeL0ref(obj,node)

    wedgeNum=size(obj.Kwedge_Vec);
    wedgeNum=wedgeNum(1);

    obj.L0ref_Mat=zeros(wedgeNum,15);

    for i=1:length(obj.Kwedge_Vec)

        node1=node.coordinates_Mat(obj.wedgeConnect_Mat(i,1),:);
        node2=node.coordinates_Mat(obj.wedgeConnect_Mat(i,2),:);
        node3=node.coordinates_Mat(obj.wedgeConnect_Mat(i,3),:);
        node4=node.coordinates_Mat(obj.wedgeConnect_Mat(i,4),:);
        node5=node.coordinates_Mat(obj.wedgeConnect_Mat(i,5),:);
        node6=node.coordinates_Mat(obj.wedgeConnect_Mat(i,6),:);

        obj.L0ref_Mat(i,1)=norm(node1-node2);        
        obj.L0ref_Mat(i,2)=norm(node2-node3);   
        obj.L0ref_Mat(i,3)=norm(node1-node3);   
        obj.L0ref_Mat(i,4)=norm(node4-node5);   
        obj.L0ref_Mat(i,5)=norm(node5-node6);   
        obj.L0ref_Mat(i,6)=norm(node4-node6);   
        obj.L0ref_Mat(i,7)=norm(node1-node4);   
        obj.L0ref_Mat(i,8)=norm(node2-node5);   
        obj.L0ref_Mat(i,9)=norm(node3-node6);   
        obj.L0ref_Mat(i,10)=norm(node1-node5);   
        obj.L0ref_Mat(i,11)=norm(node4-node2);   
        obj.L0ref_Mat(i,12)=norm(node2-node6);   
        obj.L0ref_Mat(i,13)=norm(node5-node3);   
        obj.L0ref_Mat(i,14)=norm(node3-node4);   
        obj.L0ref_Mat(i,15)=norm(node6-node1);   
    end

end
