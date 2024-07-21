%% Calculate rotation of springs
% the function calculate the rotation angle of rotational springs when 
% given nodal coordinates for of the structure


function theta=Theta(obj,X)
       
    Xi=X(1,:);
    Xj=X(2,:);
    Xk=X(3,:);

    v1=Xi-Xj;
    v2=Xk-Xj;
    
    v1=v1/norm(v1);
    v2=v2/norm(v2);

    theta=acos(dot(v1,v2));
    
end