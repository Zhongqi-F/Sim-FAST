%% Solve for the local stiffness of rotational spring
% This function use the central difference method to calculate the local
% stiffness matrix of the rotational spring element. 

function [Klocal]=Solve_Local_Stiff(obj,X,theta0,K)

    % The flat rotational spring elemement has 3 nodes, the hessian should 
    % be a 9 by 9 matrix. Local stiffness is the numerical hessian of the 
    % potential.

    Klocal=zeros(9,9);
    
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

    % gradient
    dfdxi=1/l1/l2*(xk-xj)-1/l1^3/l2*(xi-xj)*(xk-xj)'*(xi-xj);
    dfdxk=1/l1/l2*(xi-xj)-1/l2^3/l1*(xi-xj)*(xk-xj)'*(xk-xj);
    dfdxj=-dfdxi-dfdxk;

    dfdx=[dfdxi,dfdxj,dfdxk]';
    
    % hessian
    df2dxi2=3*(xi-xj)*(xk-xj)'/l1^5/l2*(xi-xj)'*(xi-xj)-...
            (xi-xj)*(xk-xj)'/l1^3/l2*eye(3)-...
            1/l1^3/l2*((xi-xj)'*(xk-xj)+(xk-xj)'*(xi-xj));
    df2dxk2=3*(xi-xj)*(xk-xj)'/l1/l2^5*(xk-xj)'*(xk-xj)-...
            (xi-xj)*(xk-xj)'/l1/l2^5*eye(3)-...
            1/l1/l2^3*((xi-xj)'*(xk-xj)+(xk-xj)'*(xi-xj));

    df2dxidxk=1/l1/l2*eye(3)-1/l1/l2^3*(xk-xj)'*(xk-xj)-...
            1/l1^3/l2*(xi-xj)'*(xi-xj)+...
            1/l1^3/l2^3*(xi-xj)*(xk-xj)'*(xi-xj)'*(xk-xj);

    df2dxkdxi=df2dxidxk';

    df2dxidxj=-df2dxi2-df2dxidxk;
    df2dxkdxj=-df2dxk2-df2dxkdxi;

    df2dxjdxi=df2dxidxj';
    df2dxjdxk=df2dxkdxj';

    df2dxj2=-df2dxidxj-df2dxkdxj;


    df2dx2=[df2dxi2 df2dxidxj df2dxidxk;
            df2dxjdxi df2dxj2 df2dxjdxk;
            df2dxkdxi df2dxkdxj df2dxk2];
    
    % theta results
    dthetadx=-1/sqrt(1-f)*dfdx;
    
    dtheta2dx2=-1/2/((1-f)^1.5)*(dfdx*dfdx')-1/sqrt(1-f)*df2dx2;

    Klocal=K*(dthetadx*dthetadx')+K*(theta-theta0)*dtheta2dx2;

end