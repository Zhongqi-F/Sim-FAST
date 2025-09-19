function  PE = Potential(obj,theta,theta0,K)

    % This is a simple linear elastic rotational spring element 
    PE=1/2*K*(theta-theta0)^2;
        
end
