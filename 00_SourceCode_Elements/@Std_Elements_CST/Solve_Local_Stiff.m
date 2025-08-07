function [Klocal]=Solve_Local_Stiff(obj,x,X,E,v,t,A)

    % The flat rotational spring elemement has 3 nodes, the hessian should 
    % be a 9 by 9 matrix. Local stiffness is the numerical hessian of the 
    % potential.
    Klocal=zeros(9,9);
    
    % this is the original location of nodes
    X1=X(1,:);
    X2=X(2,:);
    X3=X(3,:);

    % These are the original length of the side
    L1=norm(X2-X3);
    L2=norm(X1-X3);
    L3=norm(X1-X2);

    % This is the deformed configuration
    x1=x(1,:);
    x2=x(2,:);
    x3=x(3,:);

    % These are the deformed length of the triangle
    l1=norm(x2-x3);
    l2=norm(x1-x3);
    l3=norm(x1-x2);

    % We can find the engineering strain of each side 
    epsilon1= (l1-L1)/L1;
    epsilon2= (l2-L2)/L2;
    epsilon3= (l3-L3)/L3;

    eps=[epsilon1;
         epsilon2;
         epsilon3];


    % We also want to update the alpha angles and transformation matrix
    % to ensure better convergence behavior. 
    % First we solve for the beta angles
    cos1=dot((x2-x1),(x3-x1))/l2/l3;
    cos2=dot((x1-x2),(x3-x2))/l1/l3;
    cos3=dot((x1-x3),(x2-x3))/l1/l2;

    beta1=acos(cos1);
    beta2=acos(cos2);
    beta3=acos(cos3);
    
    % Solve for the transformation matrix to convert bar strain
    % to the cst strain. 
    alpha2=beta2;
    alpha3=pi-beta3;    

    tan2=tan(alpha2);
    tan3=tan(alpha3);

    cot2=cot(alpha2);
    cot3=cot(alpha3);

    sin2=sin(alpha2);
    sin3=sin(alpha3);

    cos2=cos(alpha2);
    cos3=cos(alpha3);

    % These are the factors
    B1=(tan3-tan2)/(cot2-cot3);
    B2=-1/(cot2-cot3*sin3*cos3);
    B3=1/(cot2-cot3*sin2*cos2);

    C1=(tan3*tan3-tan2*tan2)/2/(tan2-tan3);
    C2=-1/2/cos3/cos3/(tan2-tan3);
    C3=1/2/cos2/cos2/(tan2-tan3);

    trans_mat=zeros(3,3);
    trans_mat(1,:)=[1 0 0];
    trans_mat(2,1)=B1;
    trans_mat(2,2)=B2;
    trans_mat(2,3)=B3;
    trans_mat(3,1)=C1;
    trans_mat(3,2)=C2;
    trans_mat(3,3)=C3;

    de1dx=[zeros(3,1);
           (x2-x3)';
           (x3-x2)'];
    de2dx=[(x1-x3)';
           zeros(3,1);
           (x3-x1)'];
    de3dx=[(x1-x2)';
           (x2-x1)';
           zeros(3,1);];

    dedx=[de1dx de2dx de3dx];

    d2edx2=zeros(3,9,9);
    d2edx2(1,:,:)=1/L1/l1*[ zeros(3,9);
                            zeros(3), eye(3), -eye(3);
                            zeros(3), -eye(3), eye(3)]+...
                  1/L1/l1^3*[zeros(3,1);(x2-x3)';(x3-x2)']*...
                            [zeros(3,1);(x2-x3)';(x3-x2)']';
                   
    d2edx2(2,:,:)=1/L2/l2*[eye(3), zeros(3), -eye(3); 
                           zeros(3,9); 
                           -eye(3), zeros(3), eye(3)]+...
                  1/L2/l2^3*[(x1-x3)';zeros(3,1);(x3-x1)']*...
                            [(x1-x3)';zeros(3,1);(x3-x1)']';

    d2edx2(3,:,:)=1/L1/l1*[ eye(3), -eye(3), zeros(3);
                            -eye(3), eye(3), zeros(3);
                            zeros(3,9); ]+...
                  1/L1/l1^3*[(x1-x2)';(x2-x1)';zeros(3,1)]*...
                            [(x1-x2)';(x2-x1)';zeros(3,1)]';

    EMat=E/(1-v*v)*[1 v 0;
                v 1 0;
                0 0 (1-v)/2];

    strainMat=obj.Solve_Strain(x,X);

    strainVec=[strainMat(1,1);
               strainMat(2,2);
               strainMat(1,2)+strainMat(2,1)];

    Klocal=A*t*dedx*trans_mat'*EMat*trans_mat*dedx'+...
        A*t*squeeze(tensorprod((strainVec'*EMat*trans_mat),d2edx2,2,1));
    
end