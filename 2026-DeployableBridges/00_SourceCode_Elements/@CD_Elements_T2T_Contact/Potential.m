function PE=Potential(obj,d,k,d0)
    if d > d0
        PE = 0;
    elseif d == 0
        PE = inf;
    else
        % PE = k*(log(sec(pi/2-pi/2*d/d0))-1/2*(pi/2-pi/2*d/d0)^2); 
        PE = k*(-log(d/d0)*((d-d0)^2)); 
    end
end