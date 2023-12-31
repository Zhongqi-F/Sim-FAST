
function InitializeSpr(obj,node)

    obj.currentTheta_Vec = ...
    obj.Spr_Theta(node,node.currentU_Mat);
    obj.theta_StressFree_Vec = obj.currentTheta_Vec;

    sprIJKL=obj.sprIJKL_Mat;
    numSpr=size(sprIJKL);
    numSpr=numSpr(1);

    % check the ijkl assignment to make sure they are the same direction
    for i=1:numSpr
        n1=sprIJKL(i,1);
        n2=sprIJKL(i,2);
        n3=sprIJKL(i,3);
        n4=sprIJKL(i,4);

        x1=node.coordinates_Mat(n1,:);
        x2=node.coordinates_Mat(n2,:);
        x3=node.coordinates_Mat(n3,:);
        x4=node.coordinates_Mat(n4,:);

        v1=x1-x2;
        v2=x3-x2;

        norm1=cross(v1,v2);

        if norm1(3)>0
        else
            obj.sprIJKL_Mat(i,2)=n3;
            obj.sprIJKL_Mat(i,3)=n2;
        end

    end


end

