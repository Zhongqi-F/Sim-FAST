%% Constitutive relationships for springs
% This function calculate the moment and stiffness of the rotation springs
% given the properties of the rotational springs. 


function  PE = Potential(obj,theta,theta0,K)

    % This is a simple linear elastic rotational spring element 
    PE=1/2*K*(theta-theta0)^2;
        
end
