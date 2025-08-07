function [Flocal]=Solve_Local_Force(obj,x,theta0,K)

    % The rotational spring elemement has 4 nodes, the gradient should be 
    % a 12 by 1 vector. Local force is the numerical gradient of the bar 
    % potential.    
    Flocal=zeros(12,1);

    nodei=x(1,:);
    nodej=x(2,:);
    nodek=x(3,:);
    nodel=x(4,:);
       
    rij=(nodei-nodej)';
    rkj=(nodek-nodej)';
    rkl=(nodek-nodel)';
    m=cross(rij,rkj);
    n=cross(rkj,rkl);

    theta=obj.Solve_Theta(x);

    theta1=obj.theta1;
    theta2=obj.theta2;

    if theta<theta1
        ratio=(sec(pi*(theta-theta1)/2/theta1))^2;
        newK=K+K*ratio-K;
        M=K*(theta-theta0)+(2*theta1*K/pi)*tan(pi*(theta-theta1)/2/theta1)...
            -K*theta+K*theta1;

    elseif theta>theta2
        ratio=(sec(pi*(theta-theta2)/(4*pi-2*theta2)))^2;
        newK=K+K*ratio-K;
        M=K*(theta-theta0)+(2*(2*pi-theta2)*K/pi)*tan(pi*(theta-theta2)/...
            (4*pi-2*theta2))...
            -K*theta+K*theta2;
    else   
        newK=K; 
        M=K*(theta-theta0);
    end

    parti=norm(rkj)/norm(m)/norm(m)*m;
    partl=-norm(rkj)/norm(n)/norm(n)*n;
    partj=((dot(rij,rkj)/norm(rkj)/norm(rkj))-1)*parti-dot(rkl,rkj)/norm(rkj)/norm(rkj)*partl;
    partk=((dot(rkl,rkj)/norm(rkj)/norm(rkj))-1)*partl-dot(rij,rkj)/norm(rkj)/norm(rkj)*parti;

    part=[parti;partj;partk;partl];  
    Flocal=M*part;        

end