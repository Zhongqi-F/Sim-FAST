%% Calculate the strain of bars
% The fuction calculate the strains of bars when given the deformation
% information of the structures.


function [Ex]=Bar_Strain(obj,node,U)


    nodalCoordinates=node.coordinates_Mat;
    barConnect=obj.barConnect_Mat;
    barLength=obj.L_Vec;
    
    %% Non-vectorized Version

    Ex=zeros(size(barLength));
    A=size(barLength);
    N=A(1);
    
    for i=1:N
        node1=barConnect(i,1);
        node2=barConnect(i,2);
        iden=eye(3);
        tempU=[U(node1,:)';U(node2,:)'];

        B1=1/(barLength(i)^2)*[-(nodalCoordinates(node2,:)-nodalCoordinates(node1,:)) ...
            (nodalCoordinates(node2,:)-nodalCoordinates(node1,:)) ];
        B2=1/(barLength(i)^2)*[iden -iden;-iden iden];
        Ex(i)=B1*tempU+0.5*tempU'*B2*tempU;
    end
    
    
end