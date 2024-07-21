function [Kcst]=Solve_Global_Stiff(obj,U,dedx,d2edx2,cst_strain_mat,trans_mat)

    % Load values needed for calculation
    cst_ijk=obj.cst_ijk_mat;

    % Find the number of CST elements
    Ncst=size(cst_ijk);
    Ncst=Ncst(1);

    % mechanical properties
    E_vec=obj.E_vec;
    t_vec=obj.t_vec;
    v_vec=obj.v_vec;
    
    % matrix for material properties
    E_mat=zeros(Ncst,3,3);
    E_mat(:,1,1)=E_vec./(1-v_vec.*v_vec);
    E_mat(:,2,2)=E_vec./(1-v_vec.*v_vec);

    E_mat(:,1,2)=v_vec.*E_vec./(1-v_vec.*v_vec);
    E_mat(:,2,1)=v_vec.*E_vec./(1-v_vec.*v_vec);

    E_mat(:,3,3)=E_vec./(1+v_vec)/2;

    % retrieve information about the element
    nodeNum=size(U);
    nodeNum=nodeNum(1);
 
    % The stiffness comes from two part
    % We will solve them separately
    K1cst_local=zeros(Ncst,9,9);  
    K2cst_local=zeros(Ncst,9,9);  

    % Let's first work out the K1
    % This is the geometric stiffness part
    % This part uses the d2edx2 (Ncst,3,9,9)
    d2edx2_cst=zeros(size(d2edx2));

    for i=1:9
        for j=1:9
            d2edx2_cst(:,1,i,j)=sum(squeeze(trans_mat(:,1,:))...
                .*squeeze(d2edx2(:,:,i,j)),2);
            d2edx2_cst(:,2,i,j)=sum(squeeze(trans_mat(:,2,:))...
                .*squeeze(d2edx2(:,:,i,j)),2);
            d2edx2_cst(:,3,i,j)=sum(squeeze(trans_mat(:,3,:))...
                .*squeeze(d2edx2(:,:,i,j)),2);            
        end
    end

    % The stress results
    sigma_mat=zeros(Ncst,3);
    sigma_mat(:,1)=sum(squeeze(E_mat(:,1,:)).*cst_strain_mat,2);
    sigma_mat(:,2)=sum(squeeze(E_mat(:,2,:)).*cst_strain_mat,2);
    sigma_mat(:,3)=sum(squeeze(E_mat(:,3,:)).*cst_strain_mat,2);

    % Find the derivative of dUcst/de
    dude=zeros(Ncst,3);
    dude=sigma_mat.*(obj.A_vec*ones(1,3)).*(obj.t_vec*ones(1,3));
    
    for i=1:9
        for j=1:9
            K1cst_local(:,i,j)=sum(dude(:,:).*d2edx2_cst(:,:,i,j),2);    
        end
    end


    % Next, we can start working on the other component K2
    % This is the basic stiffness matrix
    % This part uses the dedx (Ncst,3,9)

    dedx_cst=zeros(size(dedx));
    % This matrix stores the cst strain, it is computed from 
    % the bar strain using the following set of codes

    for i=1:9
        dedx_cst(:,1,i)=sum(squeeze(trans_mat(:,1,:))...
            .*squeeze(dedx(:,:,i)),2);
        dedx_cst(:,2,i)=sum(squeeze(trans_mat(:,2,:))...
            .*squeeze(dedx(:,:,i)),2);
        dedx_cst(:,3,i)=sum(squeeze(trans_mat(:,3,:))...
            .*squeeze(dedx(:,:,i)),2);
    end

    % This is the second derivative of d2Ucst/de2
    d2ude2=zeros(Ncst,3,3);
   
    d2ude2(:,1,:)=squeeze(E_mat(:,1,:)).*obj.A_vec.*obj.t_vec;
    d2ude2(:,2,:)=squeeze(E_mat(:,2,:)).*obj.A_vec.*obj.t_vec;
    d2ude2(:,3,:)=squeeze(E_mat(:,3,:)).*obj.A_vec.*obj.t_vec;

    % a temporary storage matrix
    Rhs=zeros(Ncst,3,9);
    
    for i=1:9
        Rhs(:,1,i)= sum(squeeze(d2ude2(:,1,:)).*squeeze(dedx_cst(:,:,i)),2);
        Rhs(:,2,i)= sum(squeeze(d2ude2(:,2,:)).*squeeze(dedx_cst(:,:,i)),2);
        Rhs(:,3,i)= sum(squeeze(d2ude2(:,3,:)).*squeeze(dedx_cst(:,:,i)),2);
    end

    for i=1:9
        for j=1:9
            K2cst_local(:,i,j)= sum(squeeze(dedx_cst(:,:,i))...
                .*squeeze(Rhs(:,:,j)),2);
        end
    end
    
    % The final stiffness matrix
    % The size is based on the total number of nodes
    Kcst=zeros(nodeNum*3,nodeNum*3);

    % Finally, put the local matrix into the right location 
    for i=1:Ncst

        Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)=...
            Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)+...
            squeeze(K1cst_local(i,1:3,1:3))+...
            squeeze(K2cst_local(i,1:3,1:3));

        Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)=...
            Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)+...
            squeeze(K1cst_local(i,1:3,4:6))+...
            squeeze(K2cst_local(i,1:3,4:6));

        Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)=...
            Kcst(cst_ijk(i,1)*3-2:cst_ijk(i,1)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)+...
            squeeze(K1cst_local(i,1:3,7:9))+...
            squeeze(K2cst_local(i,1:3,7:9));

        Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)=...
            Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)+...
            squeeze(K1cst_local(i,4:6,1:3))+...
            squeeze(K2cst_local(i,4:6,1:3));

        Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)=...
            Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)+...
            squeeze(K1cst_local(i,4:6,4:6))+...
            squeeze(K2cst_local(i,4:6,4:6));

        Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)=...
            Kcst(cst_ijk(i,2)*3-2:cst_ijk(i,2)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)+...
            squeeze(K1cst_local(i,4:6,7:9))+...
            squeeze(K2cst_local(i,4:6,7:9));

        Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)=...
            Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,1)*3-2:cst_ijk(i,1)*3)+...
            squeeze(K1cst_local(i,7:9,1:3))+...
            squeeze(K2cst_local(i,7:9,1:3));

        Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)=...
            Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,2)*3-2:cst_ijk(i,2)*3)+...
            squeeze(K1cst_local(i,7:9,4:6))+...
            squeeze(K2cst_local(i,7:9,4:6));

        Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)=...
            Kcst(cst_ijk(i,3)*3-2:cst_ijk(i,3)*3,...
            cst_ijk(i,3)*3-2:cst_ijk(i,3)*3)+...
            squeeze(K1cst_local(i,7:9,7:9))+...
            squeeze(K2cst_local(i,7:9,7:9));
       
    end

end