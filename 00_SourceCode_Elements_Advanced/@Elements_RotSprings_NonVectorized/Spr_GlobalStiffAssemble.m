%% Assemble global stiffness for springs
% This fucntion Calculate the global stiffness matrix from the 
% rotational spring elements


function [Kspr]=Spr_GlobalStiffAssemble(obj,node,U,Mspr,Cspr)

    A=size(U);
    NodeNum=A(1);
    Kspr=zeros(3*NodeNum,3*NodeNum);

    sprIJKL=obj.sprIJKL_Mat;
    nodalCoordinates=node.coordinates_Mat;
    

    %% This part of the code is not vectorized

    A=size(sprIJKL);
    sprNum=A(1);

    for i=1:sprNum

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

        rkj_norm=norm(rkj);
        m_square=m'*m;
        n_square=n'*n;

        parti=norm(rkj)/m_square*m;
        partl=-norm(rkj)/n_square*n;

        partj=((dot(rij,rkj)/rkj_norm/rkj_norm-1))*parti-dot(rkl,rkj)/rkj_norm/rkj_norm*partl;
        partk=((dot(rkl,rkj)/rkj_norm/rkj_norm-1))*partl-dot(rij,rkj)/rkj_norm/rkj_norm*parti;

        part=[parti;partj;partk;partl];  

        A=dot(rij,rkj)/rkj_norm/rkj_norm;
        B=dot(rkl,rkj)/rkj_norm/rkj_norm;

        partAj=1/rkj_norm/rkj_norm*((2*A-1)*rkj-rij);
        partBj=1/rkj_norm/rkj_norm*(2*B*rkj-rkl);
        partAk=1/rkj_norm/rkj_norm*(-2*A*rkj+rij);
        partBk=1/rkj_norm/rkj_norm*((1-2*B)*rkj+rkl);

        part2ii=-norm(rkj)/(m_square*m_square)*(m*(cross(rkj,m))'+cross(rkj,m)*m');
        part2ll=norm(rkj)/(n_square*n_square)*(n*(cross(rkj,n))'+cross(rkj,n)*n');
        part2ik=m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rij,m))'+cross(rij,m)*m');
        part2lj=n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkl,n))'+cross(rkl,n)*n');
        part2ij=-m*rkj'/(m_square*norm(rkj))+norm(rkj)/(m_square*m_square)*(m*(cross(rkj-rij,m))'+cross(rkj-rij,m)*m');
        part2lk=-n*rkj'/(n_square*norm(rkj))-norm(rkj)/(n_square*n_square)*(n*(cross(rkj-rkl,n))'+cross(rkj-rkl,n)*n');
        part2jj=parti*partAj'+(A-1)*part2ij-(partl*partBj'+B*part2lj);
        part2jk=parti*partAk'+(A-1)*part2ik-(partl*partBk'+B*part2lk);
        part2kk=partl*partBk'+(B-1)*part2lk-(parti*partAk'+A*part2ik);       
        part2li=zeros(3);

        part2=[part2ii part2ij part2ik part2li';
        part2ij' part2jj part2jk part2lj';
        part2ik' part2jk' part2kk part2lk';
        part2li part2lj part2lk part2ll];

        localK=Cspr(i)*(part*part')+Mspr(i)*part2;

        Kspr(index1:index1+2,index1:index1+2)=Kspr(index1:index1+2,index1:index1+2)+localK(1:3,1:3);
        Kspr(index1:index1+2,index2:index2+2)=Kspr(index1:index1+2,index2:index2+2)+localK(1:3,4:6);
        Kspr(index1:index1+2,index3:index3+2)=Kspr(index1:index1+2,index3:index3+2)+localK(1:3,7:9);
        Kspr(index1:index1+2,index4:index4+2)=Kspr(index1:index1+2,index4:index4+2)+localK(1:3,10:12);

        Kspr(index2:index2+2,index1:index1+2)=Kspr(index2:index2+2,index1:index1+2)+localK(4:6,1:3);
        Kspr(index2:index2+2,index2:index2+2)=Kspr(index2:index2+2,index2:index2+2)+localK(4:6,4:6);
        Kspr(index2:index2+2,index3:index3+2)=Kspr(index2:index2+2,index3:index3+2)+localK(4:6,7:9);
        Kspr(index2:index2+2,index4:index4+2)=Kspr(index2:index2+2,index4:index4+2)+localK(4:6,10:12);

        Kspr(index3:index3+2,index1:index1+2)=Kspr(index3:index3+2,index1:index1+2)+localK(7:9,1:3);
        Kspr(index3:index3+2,index2:index2+2)=Kspr(index3:index3+2,index2:index2+2)+localK(7:9,4:6);
        Kspr(index3:index3+2,index3:index3+2)=Kspr(index3:index3+2,index3:index3+2)+localK(7:9,7:9);
        Kspr(index3:index3+2,index4:index4+2)=Kspr(index3:index3+2,index4:index4+2)+localK(7:9,10:12);

        Kspr(index4:index4+2,index1:index1+2)=Kspr(index4:index4+2,index1:index1+2)+localK(10:12,1:3);
        Kspr(index4:index4+2,index2:index2+2)=Kspr(index4:index4+2,index2:index2+2)+localK(10:12,4:6);
        Kspr(index4:index4+2,index3:index3+2)=Kspr(index4:index4+2,index3:index3+2)+localK(10:12,7:9);
        Kspr(index4:index4+2,index4:index4+2)=Kspr(index4:index4+2,index4:index4+2)+localK(10:12,10:12);

    end
end

