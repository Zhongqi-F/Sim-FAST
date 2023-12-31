
function InitializeSpr(obj,node)


    rotSprIJKL=obj.rotSprIJKL_Mat;
    numSpr=size(rotSprIJKL);
    numSpr=numSpr(1);

    % check the ijkl assignment to make sure they are the same direction
    for i=1:numSpr
        n1=rotSprIJKL(i,1);
        n2=rotSprIJKL(i,2);
        n3=rotSprIJKL(i,3);
        n4=rotSprIJKL(i,4);

        x1=node.coordinates_Mat(n1,:);
        x2=node.coordinates_Mat(n2,:);
        x3=node.coordinates_Mat(n3,:);
        x4=node.coordinates_Mat(n4,:);

        v1=x1-x2;
        v2=x3-x2;

        norm1=cross(v1,v2);

        if norm1(3)>0
        else
            obj.rotSprIJKL_Mat(i,2)=n3;
            obj.rotSprIJKL_Mat(i,3)=n2;
        end

        X=[x1;x2;x3;x4;];
        obj.theta_Current_Vec(i)=obj.Theta(X);
        obj.theta_StressFree_Vec(i)=obj.theta_Current_Vec(i);

    end


end

