%% Find the strain matrix (2x2) for the CST element
% This function itentifies the 2x2 strain matrix using the deformed nodal
% corrdinates and the undeformed nodal coordinates

function strainMat = CalcStrain(obj,x,Xoriginal)
    
    % this is the original location of nodes
    N1=Xoriginal(1,:);
    N2=Xoriginal(2,:);
    N3=Xoriginal(3,:);

    % These are the original length of the side
    L1=norm(N2-N3);
    L2=norm(N1-N3);
    L3=norm(N1-N2);

    % This is the deformed configuration
    n1=x(1,:);
    n2=x(2,:);
    n3=x(3,:);

    % These are the deformed length of the triangle
    l1=norm(n2-n3);
    l2=norm(n1-n3);
    l3=norm(n1-n2);

    % We can find the engineering strain of each side 
    epsilon1= (l1-L1)/L1;
    epsilon2= (l2-L2)/L2;
    epsilon3= (l3-L3)/L3;

    % The next job is to convert the strain along each side of the system
    % to the strain of the CST element. This is done through a linear solve

    v1=(n2-n1)/norm(n2-n1);
    v2=(n3-n1)/norm(n3-n1);
    cosa1=dot(v1,v2); 
    sina1=sqrt(1-cosa1*cosa1);

    v1=(n1-n2)/norm(n1-n2);
    v2=(n3-n2)/norm(n3-n2);
    cosa2=dot(v1,v2); 
    sina2=sqrt(1-cosa2*cosa2);

    v1=(n1-n3)/norm(n1-n3);
    v2=(n2-n3)/norm(n2-n3);
    cosa3=dot(v1,v2);
    sina3=sqrt(1-cosa3*cosa3);

    % Assuming that cosa2 and cosa3 are not zero
    % Hopefully no angle reaches 90 during the deformation
    Mat=[sina3*sina3   -2*cosa3*sina3;
         sina2*sina2   2*cosa2*sina2;];

    Rhs=[epsilon2-cosa3*cosa3*epsilon1;
        epsilon3-cosa2*cosa2*epsilon1];

    vec=Mat\Rhs;

    strainMat=zeros(2);

    strainMat(1,1)=epsilon1;
    strainMat(2,2)=vec(1);
    strainMat(1,2)=vec(2);
    strainMat(2,1)=vec(2);


    % To avoid numerical issue, make sure that the cosine is not zero when
    % performing the linear solve. 
    % if cosa2~=0 && cosa3~=0
    %     sina3=-sina3;
    % 
    %     Mat=[sina3*sina3   2*cosa3*sina3;
    %          sina2*sina2   2*cosa2*sina2;];
    % 
    %     Rhs=[epsilon2-cosa3*cosa3*epsilon1;
    %         epsilon3-cosa2*cosa2*epsilon1];
    % 
    %     vec=Mat\Rhs;
    % 
    %     strainMat=zeros(2);
    % 
    %     strainMat(1,1)=epsilon1;
    %     strainMat(2,2)=vec(1);
    %     strainMat(1,2)=vec(2);
    %     strainMat(2,1)=vec(2);
    % 
    % elseif cosa2==0
    % 
    %     sina3=-sina3;
    % 
    %     Mat=[sina3*sina3   2*cosa3*sina3;
    %          sina1*sina1   2*cosa1*sina1;];
    % 
    %     Rhs=[epsilon1-cosa3*cosa3*epsilon2;
    %         epsilon3-cosa1*cosa1*epsilon2];
    % 
    %     vec=Mat\Rhs;
    % 
    %     strainMat=zeros(2);
    % 
    %     strainMat(1,1)=epsilon2;
    %     strainMat(2,2)=vec(1);
    %     strainMat(1,2)=vec(2);
    %     strainMat(2,1)=vec(2);
    % 
    % elseif cosa3==0
    % 
    %     sina1=-sina1;
    % 
    %     Mat=[sina2*sina2   2*cosa2*sina2;
    %          sina1*sina1   2*cosa1*sina1;];
    % 
    %     Rhs=[epsilon1-cosa2*cosa2*epsilon3;
    %         epsilon2-cosa1*cosa1*epsilon3];
    % 
    %     vec=Mat\Rhs;
    % 
    %     strainMat=zeros(2);
    % 
    %     strainMat(1,1)=epsilon3;
    %     strainMat(2,2)=vec(1);
    %     strainMat(1,2)=vec(2);
    %     strainMat(2,1)=vec(2);
    % 
    % end

end