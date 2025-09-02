%% Find the rotational angle for all theta
% This code calculate the all current theta value  

function Solve_Global_Theta(obj,node,U)

    % retrieve information about the element
    nodalCoordinates=node.coordinates_mat;
    rotSprIJK=obj.node_ijk_mat;
    rotSprK=obj.rot_spr_K_vec;

    % find the number of elements
    rotSprNum=size(rotSprK);
    rotSprNum=rotSprNum(1);

    
    % solve for the internal forces using central difference
    for i=1:rotSprNum
        
        node1=rotSprIJK(i,1);
        node2=rotSprIJK(i,2);
        node3=rotSprIJK(i,3);
        
        % The nodal cooridnates of the node after adding the deformation
        X1=nodalCoordinates(node1,:)+U(node1,:);
        X2=nodalCoordinates(node2,:)+U(node2,:);
        X3=nodalCoordinates(node3,:)+U(node3,:);

        X=[X1;X2;X3];

        theta=obj.Solve_Theta(X);
        obj.theta_current_vec(i)=theta;

    end    
end