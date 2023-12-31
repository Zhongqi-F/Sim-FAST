%% Assemble stiffness of bars
% This function assembles the global stiffness matrix of the bars


function [Kbar]=Bar_GlobalStiffAssemble(obj,node,U,Sx,C)
    
    nodalCoordinates=node.coordinates_Mat;
    barConnect=obj.barConnect_Mat;
    barLength=obj.L_Vec;
    barArea=obj.A_Vec;

    A=size(nodalCoordinates);
    N_node=A(1); % Number of nodes
    Kbar=zeros(3*N_node,3*N_node);
    
    A=size(C);
    Nbar=A(1); % Number of bars
    
    for i=1:Nbar
       NodeIndex1=barConnect(i,1);
       NodeIndex2=barConnect(i,2);
       node1=nodalCoordinates(NodeIndex1,:);
       node2=nodalCoordinates(NodeIndex2,:);

       B1=1/(barLength(i)^2)*[-(node2-node1) (node2-node1)];
       iden=eye(3);
       B2=1/(barLength(i)^2)*[iden -iden; -iden iden];

       Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];   
       Ktemp=C(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp)*(B1'+B2*Utemp)'+ ...
           Sx(i)*barArea(i)*barLength(i)*B2;

       index1=3*NodeIndex1-2;
       index2=3*NodeIndex2-2;

       Kbar(index1:index1+2,index1:index1+2)=Kbar(index1:index1+2,index1:index1+2)+Ktemp(1:3,1:3);
       Kbar(index2:index2+2,index2:index2+2)=Kbar(index2:index2+2,index2:index2+2)+Ktemp(4:6,4:6);
       Kbar(index1:index1+2,index2:index2+2)=Kbar(index1:index1+2,index2:index2+2)+Ktemp(1:3,4:6);
       Kbar(index2:index2+2,index1:index1+2)=Kbar(index2:index2+2,index1:index1+2)+Ktemp(4:6,1:3);
    end

    
end