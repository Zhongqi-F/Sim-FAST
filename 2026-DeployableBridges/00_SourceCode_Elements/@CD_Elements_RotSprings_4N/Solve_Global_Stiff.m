function [Krs]=Solve_Global_Stiff(obj,node,U)
    
    % retrieve information about the element
    nodalCoordinates=node.coordinates_mat;
    rotSprIJKL=obj.node_ijkl_mat;
    theta0=obj.theta_stress_free_vec;
    rotSprK=obj.rot_spr_K_vec;

    % find the number of elements
    rotSprNum=size(rotSprK);
    rotSprNum=rotSprNum(1);
    
    % find the number of nodal coordinates
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Krs=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
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

        % Find the local stiffness matrix
        Klocal=obj.Solve_Local_Stiff(X,theta0(i),rotSprK(i));

        nodeIndex=[node1,node2,node3,node4];

        % Finally, we just place the matrix at the right location
        for j=1:4
            for k=1:4
                Krs((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))=...
                    Krs((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))...
                    +Klocal((3*j-2):(3*j),(3*k-2):(3*k));
        
            end
        end

    end    
end