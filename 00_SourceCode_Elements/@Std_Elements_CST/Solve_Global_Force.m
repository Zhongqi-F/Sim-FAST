function [Tcst]=Solve_Global_Force(obj,node,U)

    % retrieve information about the element
    nodalCoordinates=node.coordinates_mat;
    cstIJK=obj.node_ijk_mat;

    E_vec=obj.E_vec;
    t_vec=obj.t_vec;
    v_vec=obj.v_vec;
    A_vec=obj.A_vec;

    % find the number of elements
    cstNum=size(E_vec);
    cstNum=cstNum(1);

    % find the number of nodal coordinates
    A=size(nodalCoordinates);
    nodeNum=A(1);
    Tcst=zeros(3*nodeNum,1);  
    
    % solve for the internal forces using central difference
    for i=1:cstNum
        
        node1=cstIJK(i,1);
        node2=cstIJK(i,2);
        node3=cstIJK(i,3);
        
        % The nodal cooridnates of the node before deformation
        X1=nodalCoordinates(node1,:);
        X2=nodalCoordinates(node2,:);
        X3=nodalCoordinates(node3,:);

        X=[X1;X2;X3];

        % The nodal cooridnates of the node after adding the deformation
        x1=nodalCoordinates(node1,:)+U(node1,:);
        x2=nodalCoordinates(node2,:)+U(node2,:);
        x3=nodalCoordinates(node3,:)+U(node3,:);

        x=[x1;x2;x3];


        % Find the local force vector
        Flocal=obj.Solve_Local_Force(x,X,E_vec(i),v_vec(i),t_vec(i),A_vec(i));

        % These lines put the local force vector at the right location
        % in the global force vector.
        Tcst((3*node1-2):(3*node1),:)=...
            Tcst((3*node1-2):(3*node1),:)+Flocal(1:3);
        Tcst((3*node2-2):(3*node2),:)=...
            Tcst((3*node2-2):(3*node2),:)+Flocal(4:6);
        Tcst((3*node3-2):(3*node3),:)=...
            Tcst((3*node3-2):(3*node3),:)+Flocal(7:9);

    end    
end