function [Kbar]=Solve_Global_Stiff(obj,node,U)
    
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
    Kbar=zeros(3*nodeNum,3*nodeNum);  

    % calculate the stiffness matrix
    for i=1:barNum

        % Here we find the index of the nodes
        node1=barConnect(i,1);
        node2=barConnect(i,2);
        
        % Here we find the nodal coordinates
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);

        Klocal=obj.Solve_Local_Stiff(X1,X2,L0(i),barE(i),barA(i));

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