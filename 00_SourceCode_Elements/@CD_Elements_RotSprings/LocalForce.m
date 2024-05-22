
function [Flocal]=LocalForce(obj,X,theta0,K)

    % The rotational spring elemement has 4 nodes, the gradient should be 
    % a 12 by 1 vector. Local force is the numerical gradient of the bar 
    % potential.
    
    Flocal=zeros(12,1);
    
    % We use this small step delta to calculate numerical gradient
    delta=obj.delta;

    % for each element in the matrix, we do the following calculation
    for i=1:12

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
        tempXfor=X;
        tempXfor(index2,index1)=tempXfor(index2,index1)+delta;
        
        % backward step, where the coordinates is moved in the selected
        % direction by -delta
        tempXback=X;
        tempXback(index2,index1)=tempXback(index2,index1)-delta;

        % find the change in theta angle
        thetaFor=obj.Theta(tempXfor);
        thetaBack=obj.Theta(tempXback);
    
        % This is the central deference equation
        Flocal(i)=1/2/delta*(obj.Potential(thetaFor,theta0,K)-...
            obj.Potential(thetaBack,theta0,K));

    end
end