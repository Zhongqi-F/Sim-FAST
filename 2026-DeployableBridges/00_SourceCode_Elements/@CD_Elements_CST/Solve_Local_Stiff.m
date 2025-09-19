function [Klocal]=Solve_Local_Stiff(obj,x,X,E,v,t)

    % The flat rotational spring elemement has 3 nodes, the hessian should 
    % be a 9 by 9 matrix. Local stiffness is the numerical hessian of the 
    % potential.

    Klocal=zeros(9,9);
    
    % We use this small step delta to calculate numerical Hessian
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
    
        % This is the central deference equation. Please note that we can 
        % make use of the force function for this calculation. 
        Klocal(i,:)=1/2/delta*(obj.Solve_Local_Force(tempXfor,X,E,v,t)-...
            obj.Solve_Local_Force(tempXback,X,E,v,t));

    end    
end