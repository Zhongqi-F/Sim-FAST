function [bar_strain_mat,l_Mat,x_reshape,trans_mat] ...
    = Solve_Bar_Strain(obj,U,X0)
    
    x=U+X0; % Deformed node
    ijk_mat=obj.cst_ijk_mat; % Connectivity matrix

    % Number of cst elements
    Ncst=size(ijk_mat);
    Ncst=Ncst(1);

    x_reshape=zeros(Ncst,9);
    x_reshape(:,1:3)=x(ijk_mat(:,1),:);
    x_reshape(:,4:6)=x(ijk_mat(:,2),:);
    x_reshape(:,7:9)=x(ijk_mat(:,3),:);

    L1_Vec=obj.L_mat(:,1);
    L2_Vec=obj.L_mat(:,2);
    L3_Vec=obj.L_mat(:,3);

    l1_Vec=x_reshape(:,4:6)-x_reshape(:,7:9);
    l2_Vec=x_reshape(:,1:3)-x_reshape(:,7:9);
    l3_Vec=x_reshape(:,4:6)-x_reshape(:,1:3);

    l1_Vec=vecnorm(l1_Vec,2,2);
    l2_Vec=vecnorm(l2_Vec,2,2);
    l3_Vec=vecnorm(l3_Vec,2,2);

    % We can find the engineering strain of each side 
    epsilon1_Vec = (l1_Vec-L1_Vec)./L1_Vec;
    epsilon2_Vec = (l2_Vec-L2_Vec)./L2_Vec;
    epsilon3_Vec = (l3_Vec-L3_Vec)./L3_Vec;

    % Output results
    bar_strain_mat=[epsilon1_Vec,epsilon2_Vec,epsilon3_Vec];
    l_Mat=[l1_Vec,L2_Vec,L3_Vec];

    
    

    % We also want to update the alpha angles and transformation matrix
    % to ensure better convergence behavior. 
    % First we solve for the beta angles
    cos1_vec=dot((x_reshape(:,4:6)-x_reshape(:,1:3)),...
        (x_reshape(:,7:9)-x_reshape(:,1:3)),2)./l2_Vec./l3_Vec;

    cos2_vec=dot((x_reshape(:,1:3)-x_reshape(:,4:6)),...
        (x_reshape(:,7:9)-x_reshape(:,4:6)),2)./l1_Vec./l3_Vec;

    cos3_vec=dot((x_reshape(:,1:3)-x_reshape(:,7:9)),...
        (x_reshape(:,4:6)-x_reshape(:,7:9)),2)./l1_Vec./l2_Vec;

    beta1_vec=acos(cos1_vec);
    beta2_vec=acos(cos2_vec);
    beta3_vec=acos(cos3_vec);
    
    % Solve for the transformation matrix to convert bar strain
    % to the cst strain. 
    alpha2_vec=beta2_vec;
    alpha3_vec=pi-beta3_vec;

    

    tan2_vec=tan(alpha2_vec);
    tan3_vec=tan(alpha3_vec);

    cot2_vec=cot(alpha2_vec);
    cot3_vec=cot(alpha3_vec);

    sin2_vec=sin(alpha2_vec);
    sin3_vec=sin(alpha3_vec);

    cos2_vec=cos(alpha2_vec);
    cos3_vec=cos(alpha3_vec);

    % These are the factors
    B1_vec=(tan3_vec-tan2_vec)./(cot2_vec-cot3_vec);
    B2_vec=(-1./(cot2_vec-cot3_vec))./sin3_vec./cos3_vec;
    B3_vec=(1./(cot2_vec-cot3_vec))./sin2_vec./cos2_vec;

    C1_vec=(tan3_vec.*tan3_vec-tan2_vec.*tan2_vec)/2./(tan2_vec-tan3_vec);
    C2_vec=-1/2./cos3_vec./cos3_vec./(tan2_vec-tan3_vec);
    C3_vec=1/2./cos2_vec./cos2_vec./(tan2_vec-tan3_vec);

    % A factor of 2 is included for shear
    C1_vec=C1_vec*2;
    C2_vec=C2_vec*2;
    C3_vec=C3_vec*2;

    trans_mat=zeros(Ncst,3,3);
    trans_mat(:,1,:)=ones(Ncst,1)*[1 0 0];
    trans_mat(:,2,1)=B1_vec;
    trans_mat(:,2,2)=B2_vec;
    trans_mat(:,2,3)=B3_vec;
    trans_mat(:,3,1)=C1_vec;
    trans_mat(:,3,2)=C2_vec;
    trans_mat(:,3,3)=C3_vec;



end