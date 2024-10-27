function [Kzlspr]=Solve_Global_Stiff(obj,node,U)
    
    nodalCoordinates=node.coordinates_mat;
    zlsprConnect=obj.node_ij_mat;
    zlsprK=obj.k_vec;

    A=size(nodalCoordinates);
    N_node=A(1); % Number of nodes

    Kzlspr=zeros(3*N_node,3*N_node);
    
    A=size(zlsprK);
    N_zlspr=A(1); % Number of zero length springs
    
    for i=1:N_zlspr
        index1=zlsprConnect(i,1);
        index2=zlsprConnect(i,2);

        Kzlspr(3*index1-2:3*index1,3*index1-2:3*index1)=...
            Kzlspr(3*index1-2:3*index1,3*index1-2:3*index1)+zlsprK(i)*eye(3);
        Kzlspr(3*index2-2:3*index2,3*index1-2:3*index1)=...
            Kzlspr(3*index2-2:3*index2,3*index1-2:3*index1)-zlsprK(i)*eye(3);
        Kzlspr(3*index1-2:3*index1,3*index2-2:3*index2)=...
            Kzlspr(3*index1-2:3*index1,3*index2-2:3*index2)-zlsprK(i)*eye(3);
        Kzlspr(3*index2-2:3*index2,3*index2-2:3*index2)=...
            Kzlspr(3*index2-2:3*index2,3*index2-2:3*index2)+zlsprK(i)*eye(3);

    end
   
end