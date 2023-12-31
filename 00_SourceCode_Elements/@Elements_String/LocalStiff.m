%% Local stiffness matrix of bar elements
% This function calcualates the local stiffness matrix for the bar element
% using the central difference method. The stiffness is the Hessian of the
% total potential. We can make use of the "LocalForce" function and the
% central difference principle to calculate the result. 

function Klocal=LocalStiff(obj,X1,X2,L0,E,A)

    % The bar elemements has two nodes, the hessian should be a 6 by 6
    % matrix. Local stiffness is the numerical hessian of the potential.
    Klocal=zeros(6,6);
    
    % We use this small step delta to calculate numerical Hessian
    delta=obj.delta;

    % forward step, where the coodinates is moved in x direction by delta 
    % value for node 1
    TempX1for=X1;
    TempX1for(1)=TempX1for(1)+delta;
    
    % backward step, where the coodinates is moved in x direction by -delta 
    % value for node 1
    TempX1back=X1;
    TempX1back(1)=TempX1back(1)-delta;

    % This is the central deference equation. Please note that we can make 
    % use of the force function for this calculation. 
    Klocal(1,:)=1/2/delta*(obj.LocalForce(TempX1for,X2,L0,E,A)-...
        obj.LocalForce(TempX1back,X2,L0,E,A));

    % We repeat the process for other coordinates
    TempX1for=X1;
    TempX1for(2)=TempX1for(2)+delta;
    TempX1back=X1;
    TempX1back(2)=TempX1back(2)-delta;
    Klocal(2,:)=1/2/delta*(obj.LocalForce(TempX1for,X2,L0,E,A)-...
        obj.LocalForce(TempX1back,X2,L0,E,A));

    TempX1for=X1;
    TempX1for(3)=TempX1for(3)+delta;
    TempX1back=X1;
    TempX1back(3)=TempX1back(3)-delta;
    Klocal(3,:)=1/2/delta*(obj.LocalForce(TempX1for,X2,L0,E,A)-...
        obj.LocalForce(TempX1back,X2,L0,E,A));

    % for the second node
    TempX2for=X2;
    TempX2for(1)=TempX2for(1)+delta;
    TempX2back=X2;
    TempX2back(1)=TempX2back(1)-delta;
    Klocal(4,:)=1/2/delta*(obj.LocalForce(X1,TempX2for,L0,E,A)-...
        obj.LocalForce(X1,TempX2back,L0,E,A));

    TempX2for=X2;
    TempX2for(2)=TempX2for(2)+delta;
    TempX2back=X2;
    TempX2back(2)=TempX2back(2)-delta;
    Klocal(5,:)=1/2/delta*(obj.LocalForce(X1,TempX2for,L0,E,A)-...
        obj.LocalForce(X1,TempX2back,L0,E,A));

    TempX2for=X2;
    TempX2for(3)=TempX2for(3)+delta;
    TempX2back=X2;
    TempX2back(3)=TempX2back(3)-delta;
    Klocal(6,:)=1/2/delta*(obj.LocalForce(X1,TempX2for,L0,E,A)-...
        obj.LocalForce(X1,TempX2back,L0,E,A));
    
end