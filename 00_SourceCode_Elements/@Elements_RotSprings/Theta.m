%% Calculate rotation of springs
% the function calculate the rotation angle of rotational springs when 
% given nodal coordinates for of the structure


function theta=Theta(obj,X)
       
    Xi=X(1,:);
    Xj=X(2,:);
    Xk=X(3,:);
    Xl=X(4,:);

    rij=Xi-Xj;
    rkj=Xk-Xj;
    rkl=Xk-Xl;

    m=cross(rij,rkj);
    n=cross(rkj,rkl);

    yita=1;                        
    if dot(m,rkl)==0
        yita=1;
    else
        yita=sign(dot(m,rkl));
    end

    theta=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi);
    
end