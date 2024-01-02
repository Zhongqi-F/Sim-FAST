%% Assemble stiffness matrix of rotational springs
% This function assembles the global stiffness matrix of the rotational
% spring elements.

function [Krs]=GlobalStiff(obj,node,U)
    
    % retrieve information about the element
    nodalCoordinates=node.coordinates_Mat;
    rotSprIJK=obj.rotSprIJK_Mat;
    theta0=obj.theta_StressFree_Vec;
    rotSprK=obj.rotSprK_Vec;

    % find the number of elements
    rotSprNum=size(rotSprK);
    rotSprNum=rotSprNum(1);
    
    % find the number of nodal coordinates
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Krs=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
    for i=1:rotSprNum

        node1=rotSprIJK(i,1);
        node2=rotSprIJK(i,2);
        node3=rotSprIJK(i,3);
        
        % The nodal cooridnates of the node after adding the deformation
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);
        X3=nodalCoordinates(node3,:)+U(node3,:);

        X=[X1;X2;X3];

        % Find the local stiffness matrix
        Klocal=obj.LocalStiff(X,theta0(i),rotSprK(i));

        nodeIndex=[node1,node2,node3];

        % Finally, we just place the matrix at the right location
        for j=1:3
            for k=1:3
                Krs((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))=...
                    Krs((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))...
                    +Klocal((3*j-2):(3*j),(3*k-2):(3*k));
        
            end
        end

    end    
end