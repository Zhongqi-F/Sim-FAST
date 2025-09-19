function [Flocal]=Solve_Local_Force(obj,x,X,E,v,t)

    % The flat rotational spring elemement has 3 nodes, the gradient 
    % should be a 9 by 1 vector. Local force is the numerical gradient of 
    % the bar potential.
    
    Flocal=zeros(9,1);
    
    % We use this small step delta to calculate numerical gradient
    delta=obj.delta;

    % for each element in the matrix, we do the following calculation
    for i=1:9

        % find the index
        if mod(i,3)==0
            index1=3;
            index2=i/3;
        else
            index1=mod(i,3);
            index2=floor(i/3)+1;
        end
        
        % forward step, where the coordinate is moved in the selected
        % direction by +delta
        tempXfor=x;
        tempXfor(index2,index1)=tempXfor(index2,index1)+delta;
        
        % backward step, where the coordinates is moved in the selected
        % direction by -delta
        tempXback=x;
        tempXback(index2,index1)=tempXback(index2,index1)-delta;

        % find the change in theta angle
        strainFor=obj.Solve_Strain(tempXfor,X);
        strainBack=obj.Solve_Strain(tempXback,X);

        % Find the area of the triangle
        N1=X(1,:);
        N2=X(2,:);
        N3=X(3,:);

        v1=(N2-N1);
        v2=(N3-N1);

        area=norm(cross(v1,v2))/2;
                    
        % This is the central deference equation
        Flocal(i)=1/2/delta*(obj.Potential(strainFor,E,v,area,t)-...
            obj.Potential(strainBack,E,v,area,t));

    end
end