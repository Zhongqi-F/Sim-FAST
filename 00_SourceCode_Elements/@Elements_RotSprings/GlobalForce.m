%% Assemble inner force of rotatinal springs
% This code calculate the global force vector of the rotational springs
% using the central difference method. 

function [Trs]=GlobalForce(obj,node,U)

    % retrieve information about the element
    nodalCoordinates=node.coordinates_Mat;
    rotSprIJKL=obj.rotSprIJKL_Mat;
    theta0=obj.theta_StressFree_Vec;
    rotSprK=obj.rotSprK_Vec;

    % find the number of elements
    rotSprNum=size(rotSprK);
    rotSprNum=rotSprNum(1);

    % find the number of nodal coordinates
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Trs=zeros(3*nodeNum,1);  
    
    % solve for the internal forces using central difference
    for i=1:rotSprNum
        
        node1=rotSprIJKL(i,1);
        node2=rotSprIJKL(i,2);
        node3=rotSprIJKL(i,3);
        node4=rotSprIJKL(i,4);
        
        % The nodal cooridnates of the node after adding the deformation
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);
        X3=nodalCoordinates(node3,:)+U(node3,:);
        X4=nodalCoordinates(node4,:)+U(node4,:);

        X=[X1;X2;X3;X4;];

        % Find the local force vector
        Flocal=obj.LocalForce(X,theta0(i),rotSprK(i));

        % These lines put the local force vector at the right location
        % in the global force vector.
        Trs((3*node1-2):(3*node1),:)=...
            Trs((3*node1-2):(3*node1),:)+Flocal(1:3);
        Trs((3*node2-2):(3*node2),:)=...
            Trs((3*node2-2):(3*node2),:)+Flocal(4:6);
        Trs((3*node3-2):(3*node3),:)=...
            Trs((3*node3-2):(3*node3),:)+Flocal(7:9);
        Trs((3*node4-2):(3*node4),:)=...
            Trs((3*node4-2):(3*node4),:)+Flocal(10:12);

    end    
end