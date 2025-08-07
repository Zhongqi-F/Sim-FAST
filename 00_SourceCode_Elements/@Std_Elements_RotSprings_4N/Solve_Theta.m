function theta=Solve_Theta(obj,X)
       
    Xi=X(1,:);
    Xj=X(2,:);
    Xk=X(3,:);
    Xl=X(4,:);

    rij=Xi-Xj;
    rkj=Xk-Xj;
    rkl=Xk-Xl;

    m=cross(rij,rkj);
    n=cross(rkj,rkl);
                   
    if dot(m,rkl)==0
        yita=1;
    else
        yita=sign(dot(m,rkl));
    end

    theta=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi);
    
end