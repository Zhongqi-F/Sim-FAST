%% Assemble global stiffness matrix of bars
% This function assembles the global stiffness matrix of the bars. After
% solving the local stiffness matrix, we still need to assemble them to
% form a global stiffness matrix. Bascially, it is to put the value from
% the local stiffness matrix at the right location in the global matrix. 

function [Kbar]=GlobalStiff(obj,node,U)
    
    % retrieve information about the elements
    nodalCoordinates=node.coordinates_Mat;
    barConnect=obj.barConnect_Mat;
    L0=obj.L0_Vec;
    barA=obj.A_Vec;
    barE=obj.E_Vec;

    barNum=size(barA);
    barNum=barNum(1);
    
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Kbar=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
    for i=1:barNum

        % Here we find the index of the nodes
        node1=barConnect(i,1);
        node2=barConnect(i,2);
        
        % Here we find the nodal coordinates
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);

        Klocal=obj.LocalStiff(X1,X2,L0(i),barE(i),barA(i));

        % Finally, we just place the matrix at the right location
        Kbar((3*node1-2):(3*node1),(3*node1-2):(3*node1))=...
            Kbar((3*node1-2):(3*node1),(3*node1-2):(3*node1))+Klocal(1:3,1:3);
        Kbar((3*node1-2):(3*node1),(3*node2-2):(3*node2))=...
            Kbar((3*node1-2):(3*node1),(3*node2-2):(3*node2))+Klocal(1:3,4:6);
        Kbar((3*node2-2):(3*node2),(3*node1-2):(3*node1))=...
            Kbar((3*node2-2):(3*node2),(3*node1-2):(3*node1))+Klocal(4:6,1:3);
        Kbar((3*node2-2):(3*node2),(3*node2-2):(3*node2))=...
            Kbar((3*node2-2):(3*node2),(3*node2-2):(3*node2))+Klocal(4:6,4:6);

    end    
end