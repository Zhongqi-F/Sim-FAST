function [Klocal]=Solve_Local_Stiff(obj,x,theta0,K)

    % The rotational spring elemement has 4 nodes, the hessian should be 
    % a 12 by 12 matrix. Local stiffness is the numerical hessian of the 
    % potential.

    Klocal=zeros(12,12);
    
    nodei=x(1,:);
    nodej=x(2,:);
    nodek=x(3,:);
    nodel=x(4,:);

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
        
    rij=(nodei-nodej)';
    rkj=(nodek-nodej)';
    rkl=(nodek-nodel)';
    m=cross(rij,rkj);
    n=cross(rkj,rkl);

    parti=norm(rkj)/norm(m)/norm(m)*m;
    partl=-norm(rkj)/norm(n)/norm(n)*n;
    partj=((dot(rij,rkj)/norm(rkj)/norm(rkj)-1))*parti-dot(rkl,rkj)/norm(rkj)/norm(rkj)*partl;
    partk=((dot(rkl,rkj)/norm(rkj)/norm(rkj)-1))*partl-dot(rij,rkj)/norm(rkj)/norm(rkj)*parti;

    part=[parti;partj;partk;partl];  
        
    A=dot(rij,rkj)/(norm(rkj)^2);
    B=dot(rkl,rkj)/(norm(rkj)^2);
    partAj=1/(norm(rkj)^2)*((2*A-1)*rkj-rij);
    partBj=1/(norm(rkj)^2)*(2*B*rkj-rkl);
    partAk=1/(norm(rkj)^2)*(-2*A*rkj+rij);
    partBk=1/(norm(rkj)^2)*((1-2*B)*rkj+rkl);        

    part2ii=-norm(rkj)/(norm(m)^4)*(m*(cross(rkj,m))'+cross(rkj,m)*m');
    part2ll=norm(rkj)/(norm(n)^4)*(n*(cross(rkj,n))'+cross(rkj,n)*n');
    part2ik=m*rkj'/((norm(m)^2)*norm(rkj))+norm(rkj)/(norm(m)^4)*(m*(cross(rij,m))'+cross(rij,m)*m');
    part2lj=n*rkj'/((norm(n)^2)*norm(rkj))-norm(rkj)/(norm(n)^4)*(n*(cross(rkl,n))'+cross(rkl,n)*n');
    part2ij=-m*rkj'/((norm(m)^2)*norm(rkj))+norm(rkj)/(norm(m)^4)*(m*(cross(rkj-rij,m))'+cross(rkj-rij,m)*m');
    part2lk=-n*rkj'/((norm(n)^2)*norm(rkj))-norm(rkj)/(norm(n)^4)*(n*(cross(rkj-rkl,n))'+cross(rkj-rkl,n)*n');
    part2jj=parti*partAj'+(A-1)*part2ij-(partl*partBj'+B*part2lj);
    part2jk=parti*partAk'+(A-1)*part2ik-(partl*partBk'+B*part2lk);
    part2kk=partl*partBk'+(B-1)*part2lk-(parti*partAk'+A*part2ik);        
    part2li=zeros(3);
    
    part2=[part2ii part2ij part2ik part2li';
           part2ij' part2jj part2jk part2lj';
           part2ik' part2jk' part2kk part2lk';
           part2li part2lj part2lk part2ll];
    
    Klocal=newK*(part*part')+M*part2;        

end