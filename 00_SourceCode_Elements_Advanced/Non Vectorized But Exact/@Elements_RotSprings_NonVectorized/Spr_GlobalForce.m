%% Assemble inner force vector for springs
% The output Tspr is the global inner force vector for rotational springs.


function [Tspr]=Spr_GlobalForce(obj,node,U,M)
    
    A=size(U);
    NodeNum=A(1);    
    Tspr=zeros(3*NodeNum,1);

    nodalCoordinates=node.coordinates_Mat;
    sprIJKL=obj.sprIJKL_Mat;
    sprRotK=obj.sprRotK_Vec;

    
    %% The non-vectorized version of the code
    
    A=size(sprRotK);    
    BarNum=A(1);
    for i=1:BarNum

        if sprIJKL(i,1)==0   
        else
            nodei=nodalCoordinates(sprIJKL(i,1),:)+U(sprIJKL(i,1),:);
            nodej=nodalCoordinates(sprIJKL(i,2),:)+U(sprIJKL(i,2),:);
            nodek=nodalCoordinates(sprIJKL(i,3),:)+U(sprIJKL(i,3),:);
            nodel=nodalCoordinates(sprIJKL(i,4),:)+U(sprIJKL(i,4),:);

            index1=3*sprIJKL(i,1)-2;
            index2=3*sprIJKL(i,2)-2;
            index3=3*sprIJKL(i,3)-2;
            index4=3*sprIJKL(i,4)-2;

            rij=(nodei-nodej)';
            rkj=(nodek-nodej)';
            rkl=(nodek-nodel)';
            m=cross(rij,rkj);
            n=cross(rkj,rkl);

            parti=norm(rkj)/norm(m)/norm(m)*m;
            partl=-norm(rkj)/norm(n)/norm(n)*n;
            partj=((dot(rij,rkj)/norm(rkj)/norm(rkj))-1)*parti-dot(rkl,rkj)/norm(rkj)/norm(rkj)*partl;
            partk=((dot(rkl,rkj)/norm(rkj)/norm(rkj))-1)*partl-dot(rij,rkj)/norm(rkj)/norm(rkj)*parti;

            part=[parti;partj;partk;partl];  
            localT=M(i)*part;    

            Tspr(index1:index1+2)=Tspr(index1:index1+2)+localT(1:3);
            Tspr(index2:index2+2)=Tspr(index2:index2+2)+localT(4:6);
            Tspr(index3:index3+2)=Tspr(index3:index3+2)+localT(7:9);
            Tspr(index4:index4+2)=Tspr(index4:index4+2)+localT(10:12);

        end
    end




end
