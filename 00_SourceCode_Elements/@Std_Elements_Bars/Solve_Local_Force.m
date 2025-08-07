function Flocal=Solve_Local_Force(obj,x1,x2,L0,E,A)

    % The bar elemements has two nodes, the gradient should be a 6 by 1
    % vector. Local force is the numerical gradient of the bar potential.
    Flocal=zeros(6,1);

    % Deformed length of the bar
    l=norm(x1-x2);

    % The force vector
    Flocal=E*A/L0*(l-L0)/l*[x1-x2,x2-x1]';
    
end