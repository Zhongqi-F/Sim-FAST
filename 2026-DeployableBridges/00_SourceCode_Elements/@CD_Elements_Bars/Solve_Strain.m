function [strain_vec]=Solve_Strain(obj,node,U) 

    % Retrieve information about the bars
    nodalCoordinates=node.coordinates_mat;
    barConnect=obj.node_ij_mat;
    L0=obj.L0_vec;

    barNum=size(L0);
    barNum=barNum(1);
    strain_vec=zeros(barNum,1);

    % Compute the strain of the element
    for i=1:barNum

        node1=barConnect(i,1);
        node2=barConnect(i,2);
        
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);

        strain_vec(i)=...
            (norm(X1-X2)-L0(i))/L0(i);

    end
end