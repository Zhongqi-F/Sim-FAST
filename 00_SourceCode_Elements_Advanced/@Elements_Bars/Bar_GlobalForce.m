%% Assemble inner force of bars
% This function assembles the global inner force vector of bars 


function [Tbar]=Bar_GlobalForce(obj,node,U,Sx)

  
    nodalCoordinates=node.coordinates_Mat;
    barConnect=obj.barConnect_Mat;
    barLength=obj.L0_Vec;
    barArea=obj.A_Vec;

    A=size(nodalCoordinates);
    N=A(1);
    Tbar=zeros(3*N,1);  

    %% Vectorized code

    NodeIndex1=barConnect(:,1);
    NodeIndex2=barConnect(:,2);
    node1=nodalCoordinates(NodeIndex1,:);
    node2=nodalCoordinates(NodeIndex2,:);
    
    B1n=(1./(barLength.*barLength)).*[-(node2-node1) (node2-node1)];  
    
    iden=eye(3);
    idenMat=[iden -iden; -iden iden];
    
    Utemp1=[U(NodeIndex1,:)';U(NodeIndex2,:)'];     
    B2Utemp=(1./(barLength.*barLength)).*(idenMat*Utemp1)';
    % Vectroized operation
    
    index1=3*NodeIndex1-2;
    index2=3*NodeIndex2-2;      

    Ttemp=Sx.*barArea.*barLength;
    Ttemp=Ttemp.*(B1n+B2Utemp);
    
    index=[index1,index1+1,index1+2,index2,index2+1,index2+2];
    index=index(:);
    Ttemp=Ttemp(:);    
    
    for i=1:length(index)
        Tbar(index(i))=Tbar(index(i))+Ttemp(i);
    end
    
end