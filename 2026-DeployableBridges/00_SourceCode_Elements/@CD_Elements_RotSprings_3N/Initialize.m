%% This funciton initialize the rotational spring elements
% The initialization process include calculating the current folding angle
% and the current stress free angle of the spring elements. We set the
% current state to be the stress-free state. 

function Initialize(obj,node)

    rotSprIJK=obj.node_ijk_mat;
    numSpr=size(rotSprIJK);
    numSpr=numSpr(1);

    obj.theta_current_vec=zeros(numSpr,1);
    obj.theta_stress_free_vec=zeros(numSpr,1);

    for i=1:numSpr

        n1=rotSprIJK(i,1);
        n2=rotSprIJK(i,2);
        n3=rotSprIJK(i,3);

        x1=node.coordinates_mat(n1,:);
        x2=node.coordinates_mat(n2,:);
        x3=node.coordinates_mat(n3,:);

        X=[x1;x2;x3];
        obj.theta_current_vec(i)=obj.Solve_Theta(X);
        obj.theta_stress_free_vec(i)=obj.theta_current_vec(i);

    end
end

