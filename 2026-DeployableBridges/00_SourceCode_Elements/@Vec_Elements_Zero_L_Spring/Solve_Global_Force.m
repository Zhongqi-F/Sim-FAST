function [Tzlspr]=Solve_Global_Force(obj,node,U)
  
    nodalCoordinates=node.coordinates_mat;
    zlsprConnect=obj.node_ij_mat;
    zlsprK=obj.k_vec;

    A=size(nodalCoordinates);
    N=A(1);
    Tzlspr=zeros(3*N,1);  

    A=size(zlsprK);
    Nzlspr=A(1);

    for i=1:Nzlspr

        index1=zlsprConnect(i,1);
        index2=zlsprConnect(i,2);

        node1=nodalCoordinates(index1,:)+U(index1,:);
        node2=nodalCoordinates(index2,:)+U(index2,:);
        
        Tzlspr(3*index1-2:3*index1)=Tzlspr(3*index1-2:3*index1)...
            +zlsprK(i)*(node1'-node2');
        Tzlspr(3*index2-2:3*index2)=Tzlspr(3*index2-2:3*index2)...
            +zlsprK(i)*(node2'-node1');
    end
    
end