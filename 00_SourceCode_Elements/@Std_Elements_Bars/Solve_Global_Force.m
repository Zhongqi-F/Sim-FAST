function [Tbar]=Solve_Global_Force(obj,node,U)

    % retrieve information about the elements
    nodalCoordinates=node.coordinates_mat;
    barConnect=obj.node_ij_mat;
    L0=obj.L0_vec;
    barA=obj.A_vec;
    barE=obj.E_vec;

    barNum=size(barA);
    barNum=barNum(1);

    A=size(nodalCoordinates);
    nodeNum=A(1);
    Tbar=zeros(3*nodeNum,1);  
    
    % solve for the internal forces using central difference
    for i=1:barNum
        
        % Here we find the index of the nodes
        node1=barConnect(i,1);
        node2=barConnect(i,2);
        
        % Here we find the nodal coordinates
        x1=nodalCoordinates(node1,:)+U(node1,:);
        x2=nodalCoordinates(node2,:)+U(node2,:);

        Flocal=obj.Solve_Local_Force(x1,x2,L0(i),barE(i),barA(i));

        % These two lines put the local force vector at the right location
        % in the global force vector.
        Tbar((3*node1-2):(3*node1),:)=...
            Tbar((3*node1-2):(3*node1),:)+Flocal(1:3);
        Tbar((3*node2-2):(3*node2),:)=...
            Tbar((3*node2-2):(3*node2),:)+Flocal(4:6);

    end    
end