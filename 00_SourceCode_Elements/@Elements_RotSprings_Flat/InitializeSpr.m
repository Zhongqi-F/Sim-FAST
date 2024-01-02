%% This funciton initialize the rotational spring elements
% The initialization process include calculating the current folding angle
% and the current stress free angle of the spring elements. We set the
% current state to be the stress-free state. 

function InitializeSpr(obj,node)

    rotSprIJK=obj.rotSprIJK_Mat;
    numSpr=size(rotSprIJK);
    numSpr=numSpr(1);

    for i=1:numSpr

        n1=rotSprIJK(i,1);
        n2=rotSprIJK(i,2);
        n3=rotSprIJK(i,3);

        x1=node.coordinates_Mat(n1,:);
        x2=node.coordinates_Mat(n2,:);
        x3=node.coordinates_Mat(n3,:);

        X=[x1;x2;x3];
        obj.theta_Current_Vec(i)=obj.Theta(X);
        obj.theta_StressFree_Vec(i)=obj.theta_Current_Vec(i);

    end
end

