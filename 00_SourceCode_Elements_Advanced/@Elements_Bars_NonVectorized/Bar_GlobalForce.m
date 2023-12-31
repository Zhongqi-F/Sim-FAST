%% Assemble inner force of bars
% This function assembles the global inner force vector of bars 


function [Tbar]=Bar_GlobalForce(obj,node,U,Sx)
  
    nodalCoordinates=node.coordinates_Mat;
    barConnect=obj.barConnect_Mat;
    barLength=obj.L_Vec;
    barArea=obj.A_Vec;

    A=size(nodalCoordinates);
    N=A(1);
    Tbar=zeros(3*N,1);  

 
    for i=1:N
       NodeIndex1=barConnect(i,1);
       NodeIndex2=barConnect(i,2);
       node1=nodalCoordinates(NodeIndex1,:);
       node2=nodalCoordinates(NodeIndex2,:);

       B1=1/(barLength(i)^2)*[-(node2-node1) (node2-node1)];
       B2=1/(barLength(i)^2)*[eye(3) -eye(3); -eye(3) eye(3)];

       Utemp=[U(NodeIndex1,:)';U(NodeIndex2,:)'];   

       Ttemp=Sx(i)*barArea(i)*barLength(i)*(B1'+B2*Utemp);
       index1=3*NodeIndex1-2;
       index2=3*NodeIndex2-2;

       Tbar(index1:index1+2)=Tbar(index1:index1+2)+Ttemp(1:3);
       Tbar(index2:index2+2)=Tbar(index2:index2+2)+Ttemp(4:6);   
    end
    
end