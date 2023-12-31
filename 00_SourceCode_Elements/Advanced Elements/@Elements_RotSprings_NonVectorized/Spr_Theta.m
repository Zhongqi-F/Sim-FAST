%% Calculate rotation of springs
% the function calculate the rotation angle of rotational springs when 
% given the deformed configuration of the structure


function [theta]=Spr_Theta(obj,node,U)

    nodalCoordinate=node.coordinates_Mat;
    sprIJKL=obj.sprIJKL_Mat;

    A=size(sprIJKL);
    N=A(1);
    theta=zeros(N,1);

    %% Non-vectorized version
%     
%     for i=1:N
%         if sprIJKL(i,1)==0
%         else        
%             nodei=newNode(sprIJKL(i,1),:)+U(sprIJKL(i,1),:);
%             nodej=newNode(sprIJKL(i,2),:)+U(sprIJKL(i,2),:);
%             nodek=newNode(sprIJKL(i,3),:)+U(sprIJKL(i,3),:);
%             nodel=newNode(sprIJKL(i,4),:)+U(sprIJKL(i,4),:);
%             rij=nodei-nodej;
%             rkj=nodek-nodej;
%             rkl=nodek-nodel;
%             m=cross(rij,rkj);
%             n=cross(rkj,rkl);
%             yita=1;                        
%             if dot(m,rkl)==0
%                 yita=1;
%             else
%                 yita=sign(dot(m,rkl));
%             end
% 
%             
%             theta(i)=mod(yita*real(acos(dot(m,n)/norm(m)/norm(n))),2*pi);
%         end
%     end
    
    %% Vectorized version   

    spr_i=sprIJKL(:,1);
    spr_j=sprIJKL(:,2);
    spr_k=sprIJKL(:,3);
    spr_l=sprIJKL(:,4);

    nonzero=find(spr_i);
    spr_Num=length(nonzero);
    
    nodei=nodalCoordinate(spr_i(nonzero),:)+U(spr_i(nonzero),:);
    nodej=nodalCoordinate(spr_j(nonzero),:)+U(spr_j(nonzero),:);
    nodek=nodalCoordinate(spr_k(nonzero),:)+U(spr_k(nonzero),:);
    nodel=nodalCoordinate(spr_l(nonzero),:)+U(spr_l(nonzero),:);
    
    rij=nodei-nodej;
    rkj=nodek-nodej;
    rkl=nodek-nodel;
    
    m=cross(rij,rkj,2);
    n=cross(rkj,rkl,2);
    
    d_m_kl=dot(m,rkl,2);
    zero_index=find(d_m_kl==0);

    yita=sign(d_m_kl);
    yita(zero_index)=1;
    
    m_norm=sqrt(dot(m,m,2));
    n_norm=sqrt(dot(n,n,2));
    
    theta_local=mod(yita.*real(acos(dot(m,n,2)./m_norm./n_norm)),2*pi);
   
    for i=1:length(theta_local)
        theta(nonzero(i))=theta_local(i);
    end
    
end