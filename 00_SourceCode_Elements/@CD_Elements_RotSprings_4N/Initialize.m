%% This funciton initialize the rotational spring elements
% The initialization process include calculating the current folding angle
% and the current stress free angle of the spring elements. We set the
% current state to be the stress-free state. 

function Initialize(obj,node)

    rotSprIJKL=obj.node_ijkl_mat;
    numSpr=size(rotSprIJKL);
    numSpr=numSpr(1);

    % check the ijkl assignment to make sure they are the same direction
    for i=1:numSpr

        n1=rotSprIJKL(i,1);
        n2=rotSprIJKL(i,2);
        n3=rotSprIJKL(i,3);
        n4=rotSprIJKL(i,4);

        x1=node.coordinates_mat(n1,:);
        x2=node.coordinates_mat(n2,:);
        x3=node.coordinates_mat(n3,:);
        x4=node.coordinates_mat(n4,:);

        v1=x1-x2;
        v2=x3-x2;

        norm1=cross(v1,v2);

        if norm1(3)>0
        else
            obj.node_ijkl_mat(i,2)=n3;
            obj.node_ijkl_mat(i,3)=n2;
        end

        X=[x1;x2;x3;x4;];
        obj.theta_current_vec(i)=obj.Solve_Theta(X);
        obj.theta_stress_free_vec(i)=obj.theta_current_vec(i);

    end
end

