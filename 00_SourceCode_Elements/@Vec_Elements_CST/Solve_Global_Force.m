function [Tcst]=Solve_Global_Force(obj,U,dedx,cst_strain_mat,trans_mat)

    cst_ijf=obj.cst_ijk_mat;

    Ncst=size(cst_ijf);
    Ncst=Ncst(1);

    % mechanical properties
    E_vec=obj.E_vec;
    v_vec=obj.v_vec;

    E_mat=zeros(Ncst,3,3);
    E_mat(:,1,1)=E_vec./(1-v_vec.*v_vec);
    E_mat(:,2,2)=E_vec./(1-v_vec.*v_vec);

    E_mat(:,1,2)=v_vec.*E_vec./(1-v_vec.*v_vec);
    E_mat(:,2,1)=v_vec.*E_vec./(1-v_vec.*v_vec);

    E_mat(:,3,3)=E_vec./(1+v_vec)/2;

    % Calculate the stress based on the strain
    sigma_mat=zeros(Ncst,3);
    sigma_mat(:,1)=sum(squeeze(E_mat(:,1,:)).*cst_strain_mat,2);
    sigma_mat(:,2)=sum(squeeze(E_mat(:,2,:)).*cst_strain_mat,2);
    sigma_mat(:,3)=sum(squeeze(E_mat(:,3,:)).*cst_strain_mat,2);

    % Modify the dedx matrix with the transformation
    depdx=zeros(size(dedx));

    % first row
    depdx(:,1,1)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,1)),2);
    depdx(:,1,2)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,2)),2);
    depdx(:,1,3)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,3)),2);

    depdx(:,1,4)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,4)),2);
    depdx(:,1,5)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,5)),2);
    depdx(:,1,6)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,6)),2);

    depdx(:,1,7)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,7)),2);
    depdx(:,1,8)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,8)),2);
    depdx(:,1,9)=sum(squeeze(trans_mat(:,1,:)).*squeeze(dedx(:,:,9)),2);

    % second row
    depdx(:,2,1)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,1)),2);
    depdx(:,2,2)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,2)),2);
    depdx(:,2,3)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,3)),2);

    depdx(:,2,4)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,4)),2);
    depdx(:,2,5)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,5)),2);
    depdx(:,2,6)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,6)),2);

    depdx(:,2,7)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,7)),2);
    depdx(:,2,8)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,8)),2);
    depdx(:,2,9)=sum(squeeze(trans_mat(:,2,:)).*squeeze(dedx(:,:,9)),2);

    % third row
    depdx(:,3,1)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,1)),2);
    depdx(:,3,2)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,2)),2);
    depdx(:,3,3)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,3)),2);

    depdx(:,3,4)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,4)),2);
    depdx(:,3,5)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,5)),2);
    depdx(:,3,6)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,6)),2);

    depdx(:,3,7)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,7)),2);
    depdx(:,3,8)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,8)),2);
    depdx(:,3,9)=sum(squeeze(trans_mat(:,3,:)).*squeeze(dedx(:,:,9)),2);

    % Solve for the global force
    localF_mat=zeros(Ncst,9);

    localF_mat(:,1)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,1)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,2)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,2)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,3)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,3)),2)...
        .*obj.A_vec.*obj.t_vec;

    localF_mat(:,4)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,4)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,5)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,5)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,6)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,6)),2)...
        .*obj.A_vec.*obj.t_vec;

    localF_mat(:,7)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,7)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,8)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,8)),2)...
        .*obj.A_vec.*obj.t_vec;
    localF_mat(:,9)=sum(sigma_mat(:,:).*squeeze(depdx(:,:,9)),2)...
        .*obj.A_vec.*obj.t_vec;
    
    % find the number of nodal coordinates
    A=size(U);
    nodeNum=A(1);
    Tcst=zeros(1,3*nodeNum);  

    % Put local force values to the global matrix
    ijk_mat=obj.cst_ijk_mat;
    ijk_mat=reshape(ijk_mat',1,Ncst*3);
    localF_mat=reshape(localF_mat',1,Ncst*9);

    for i=1:Ncst*3
        nodeNum=ijk_mat(i);
        Tcst((3*nodeNum-2):3*nodeNum)=...
            Tcst((3*nodeNum-2):3*nodeNum)+localF_mat((3*i-2):3*i);
    end

    Tcst=Tcst';
end