%% Assemble inner force of wedges
% This function assembles the global inner force vector of wedges 


function [Tw]=Wedge_GlobalForce(obj,node,U)

    nodalCoordinates=node.coordinates_Mat;
    wedgeConnect=obj.wedgeConnect_Mat;
    L0ref=obj.L0ref_Mat;
    Kwedge=obj.Kwedge_Vec;

    wedgeNum=size(Kwedge);
    wedgeNum=wedgeNum(1);

    A=size(nodalCoordinates);
    nodeNum=A(1);
    Tw=zeros(3*nodeNum,1);  
    
    % solve for the internal forces using central difference
    for i=1:wedgeNum
        
        node1=wedgeConnect(i,1);
        node2=wedgeConnect(i,2);
        node3=wedgeConnect(i,3);
        node4=wedgeConnect(i,4);
        node5=wedgeConnect(i,5);
        node6=wedgeConnect(i,6);
        
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);
        X3=nodalCoordinates(node3,:)+U(node3,:);
        X4=nodalCoordinates(node4,:)+U(node4,:);
        X5=nodalCoordinates(node5,:)+U(node5,:);
        X6=nodalCoordinates(node6,:)+U(node6,:);

        X=[X1;X2;X3;X4;X5;X6];

        Flocal=obj.Wedge_LocalForce(X,Kwedge(i),L0ref(i,:));

        % These two lines put the local force vector at the right location
        % in the global force vector.
        Tw((3*node1-2):(3*node1),:)=...
            Tw((3*node1-2):(3*node1),:)+Flocal(1:3);
        Tw((3*node2-2):(3*node2),:)=...
            Tw((3*node2-2):(3*node2),:)+Flocal(4:6);
        Tw((3*node3-2):(3*node3),:)=...
            Tw((3*node3-2):(3*node3),:)+Flocal(7:9);
        Tw((3*node4-2):(3*node4),:)=...
            Tw((3*node4-2):(3*node4),:)+Flocal(10:12);
        Tw((3*node5-2):(3*node5),:)=...
            Tw((3*node5-2):(3*node5),:)+Flocal(13:15);
        Tw((3*node6-2):(3*node6),:)=...
            Tw((3*node6-2):(3*node6),:)+Flocal(16:18);

    end    
end