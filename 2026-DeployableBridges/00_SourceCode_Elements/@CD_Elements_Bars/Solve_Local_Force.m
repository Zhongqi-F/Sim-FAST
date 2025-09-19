function Flocal=Solve_Local_Force(obj,X1,X2,L0,E,A)

    % The bar elemements has two nodes, the gradient should be a 6 by 1
    % vector. Local force is the numerical gradient of the bar potential.
    Flocal=zeros(6,1);
    
    % We use this small step delta to calculate numerical gradient
    delta=obj.delta;

    % forward step, where the coodinates is moved in x direction by delta 
    % value for node 1
    TempX1for=X1;
    TempX1for(1)=TempX1for(1)+delta;
    
    % backward step, where the coodinates is moved in x direction by -delta 
    % value for node 1
    TempX1back=X1;
    TempX1back(1)=TempX1back(1)-delta;

    % This is the central deference equation
    Flocal(1)=1/2/delta*(obj.Potential(TempX1for,X2,L0,E,A)-...
        obj.Potential(TempX1back,X2,L0,E,A));

    % We repeat the process for other coordinates
    TempX1for=X1;
    TempX1for(2)=TempX1for(2)+delta;
    TempX1back=X1;
    TempX1back(2)=TempX1back(2)-delta;
    Flocal(2)=1/2/delta*(obj.Potential(TempX1for,X2,L0,E,A)-...
        obj.Potential(TempX1back,X2,L0,E,A));

    TempX1for=X1;
    TempX1for(3)=TempX1for(3)+delta;
    TempX1back=X1;
    TempX1back(3)=TempX1back(3)-delta;
    Flocal(3)=1/2/delta*(obj.Potential(TempX1for,X2,L0,E,A)-...
        obj.Potential(TempX1back,X2,L0,E,A));

    % for the second node
    TempX2for=X2;
    TempX2for(1)=TempX2for(1)+delta;
    TempX2back=X2;
    TempX2back(1)=TempX2back(1)-delta;
    Flocal(4)=1/2/delta*(obj.Potential(X1,TempX2for,L0,E,A)-...
        obj.Potential(X1,TempX2back,L0,E,A));

    TempX2for=X2;
    TempX2for(2)=TempX2for(2)+delta;
    TempX2back=X2;
    TempX2back(2)=TempX2back(2)-delta;
    Flocal(5)=1/2/delta*(obj.Potential(X1,TempX2for,L0,E,A)-...
        obj.Potential(X1,TempX2back,L0,E,A));

    TempX2for=X2;
    TempX2for(3)=TempX2for(3)+delta;
    TempX2back=X2;
    TempX2back(3)=TempX2back(3)-delta;
    Flocal(6)=1/2/delta*(obj.Potential(X1,TempX2for,L0,E,A)-...
        obj.Potential(X1,TempX2back,L0,E,A));
end