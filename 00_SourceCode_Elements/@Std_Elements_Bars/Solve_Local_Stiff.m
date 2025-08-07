function Klocal=Solve_Local_Stiff(obj,x1,x2,L0,E,A)

    % The bar elemements has two nodes, the hessian should be a 6 by 6
    % matrix. Local stiffness is the numerical hessian of the potential.
    Klocal=zeros(6,6);
    
    % The deformed length of bars
    l=norm(x1-x2);

    % Compute stiffness
    K1=E*A/l*(1/l/l)*[x1-x2,x2-x1]'*[x1-x2,x2-x1];
    K2=E*A/l*(l-L0)/L0*[eye(3) -eye(3); -eye(3) eye(3)];

    Klocal=K1+K2;
    
end