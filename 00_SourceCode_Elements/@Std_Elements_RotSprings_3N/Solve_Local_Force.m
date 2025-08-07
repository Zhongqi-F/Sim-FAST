
function [Flocal]=Solve_Local_Force(obj,X,theta0,K)

    % The flat rotational spring elemement has 3 nodes, the gradient 
    % should be a 9 by 1 vector. Local force is the numerical gradient of 
    % the bar potential.    
    Flocal=zeros(9,1);
    
    % Get the nodal coordinates
    xi=X(1,:);
    xj=X(2,:);
    xk=X(3,:);

    % Obtain the angle
    theta=obj.Solve_Theta(X);

    % Length of the two edges
    l1=sqrt((xi-xj)*(xi-xj)');
    l2=sqrt((xk-xj)*(xk-xj)');

    % Value of f and its gradient and hessian
    f=(xi-xj)*(xk-xj)'/l1/l2;

    dfdxi=1/l1/l2*(xk-xj)-1/l1^3/l2*(xi-xj)*(xk-xj)'*(xi-xj);
    dfdxk=1/l1/l2*(xi-xj)-1/l2^3/l1*(xi-xj)*(xk-xj)'*(xk-xj);
    dfdxj=-dfdxi-dfdxk;

    dfdx=[dfdxi,dfdxj,dfdxk]';

    % value of dthetadx
    dthetadx=-1/sqrt(1-f)*dfdx;

    % value of Flocal
    Flocal=K*(theta-theta0)*dthetadx; 

end