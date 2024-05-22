%% Assemble stiffness of bars
% This function assembles the global stiffness matrix of the bars

function [Kw]=Wedge_GlobalStiff(obj,node,U)
    
    nodalCoordinates=node.coordinates_Mat;
    wedgeConnect=obj.wedgeConnect_Mat;
    L0ref=obj.L0ref_Mat;
    E=obj.E_Vec;
    A0ref=obj.A0ref_Vec;

    wedgeNum=size(E);
    wedgeNum=wedgeNum(1);
    
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Kw=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
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

        Klocal=obj.Wedge_LocalStiff(X,E(i),A0ref(i),L0ref(i,:));

        nodeIndex=[node1,node2,node3,node4,node5,node6];

        % Finally, we just place the matrix at the right location
        for j=1:6
            for k=1:6
                Kw((3*nodeIndex(j)-2):(3*nodeIndex(j)),(3*nodeIndex(k)-2):(3*nodeIndex(k)))=...
                    Kw((3*nodeIndex(j)-2):(3*nodeIndex(j)),(3*nodeIndex(k)-2):(3*nodeIndex(k)))...
                    +Klocal((3*j-2):(3*j),(3*k-2):(3*k));
        
            end
        end

    end    
end