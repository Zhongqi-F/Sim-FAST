function [Kcst]=Solve_Global_Stiff(obj,node,U)
    
    % retrieve information about the element
    nodalCoordinates=node.coordinates_mat;
    cstIJK=obj.node_ijk_mat;

    E_vec=obj.E_vec;
    t_vec=obj.t_vec;
    v_vec=obj.v_vec;


    % find the number of elements
    rotSprNum=size(E_vec);
    rotSprNum=rotSprNum(1);
    
    % find the number of nodal coordinates
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Kcst=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
    for i=1:rotSprNum

        node1=cstIJK(i,1);
        node2=cstIJK(i,2);
        node3=cstIJK(i,3);
        
        % The nodal cooridnates of the node before deformation
        X1=nodalCoordinates(node1,:);
        X2=nodalCoordinates(node2,:);
        X3=nodalCoordinates(node3,:);

        Xoriginal=[X1;X2;X3];

        % The nodal cooridnates of the node after adding the deformation
        x1=nodalCoordinates(node1,:)+U(node1,:);
        x2=nodalCoordinates(node2,:)+U(node2,:);
        x3=nodalCoordinates(node3,:)+U(node3,:);

        x=[x1;x2;x3];

        % Find the local stiffness matrix
        Klocal=obj.Solve_Local_Stiff(x,Xoriginal,E_vec(i),v_vec(i),t_vec(i));

        nodeIndex=[node1,node2,node3];

        % Finally, we just place the matrix at the right location
        for j=1:3
            for k=1:3
                Kcst((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))=...
                    Kcst((3*nodeIndex(j)-2):(3*nodeIndex(j)),...
                    (3*nodeIndex(k)-2):(3*nodeIndex(k)))...
                    +Klocal((3*j-2):(3*j),(3*k-2):(3*k));
        
            end
        end

    end    
end