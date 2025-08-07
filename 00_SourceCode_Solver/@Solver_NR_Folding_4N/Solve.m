%% Nonlinear solver for folding (Assembling)

function [Uhis]=Solve(obj)
    
    increStep=obj.increStep;
    tol=obj.tol;
    iterMax=obj.iterMax;    
    supp=obj.supp;  

    assembly=obj.assembly;
    U=assembly.node.current_U_mat;

    A=size(U);
    NodeNum=A(1);
    Uhis=zeros(increStep,NodeNum,3);

    % Find the external forces that is currently applied on the structure
    currentAppliedForce=zeros(3*NodeNum,1);    
    for i=1:NodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = assembly.node.current_ext_force_mat(i,:);
    end  


   
    % Analysis start
    fprintf('Self Assemble Analysis Start');
    count=1;

    % find the zero strain angle, before and after the analysis
    sprZeroStrain_before=assembly.rot_spr_4N.theta_current_vec;
    sprZeroStrain_after=obj.targetRot;
    
        
    for i=1:increStep
        
        sprZeroStrain_current=i/increStep*sprZeroStrain_after+...
            (1-i/increStep)*sprZeroStrain_before;

        fprintf('Icrement = %d\n',count);
        step=1;
        R=1;
             
        while and(step<iterMax,R>tol)
            
            assembly.rot_spr_4N.theta_stress_free_vec=...
                sprZeroStrain_current;
            [T,K]=assembly.Solve_FK(U);

            % calculate the unbalanced force
            unbalance=currentAppliedForce-T; 
            
            [K,unbalance]=Mod_K_For_Supp(K,supp,unbalance);
            K=sparse(K);          

            deltaU=(unbalance'/K)';
            for j=1:NodeNum
                U((j),:)=U((j),:)+deltaU(3*j-2:3*j)';
            end 

            R=norm(unbalance);
            fprintf('	Iteration = %d, R = %e \n',step,R);
            step=step+1; 

        end
        
        Uhis(i,:,:)=U;
        count=count+1;
    end

    assembly.node.current_U_mat=U;
end

