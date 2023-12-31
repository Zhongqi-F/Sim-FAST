%% Potential energy of bar elements
% This function defines the potential energy of the bar element. We use a
% linear elastic (large deformation) formulation for the bar, where we
% assume that the total potential is quadratic to the extention of the bar.

function PE=Potential(obj,X1,X2,L0,E,A)

    if norm(X1-X2) > L0
        PE=0;
    else
        PE=1/2*E*A/L0*(norm(X1-X2)-L0)^2;    
    end
    % The stiffness factor is EA/L0
    % the extension is norm(X1-X2)-L0
    
end