%% Newton-Raphson solver for support deformation

function [Uhis]=Solve(obj)

    % initialize and set up storage matrix for Uhis
    increStep=obj.increStep;
    tol=obj.tol;
    iterMax=obj.iterMax;
    
    supp=obj.supp;
    suppTarget=obj.suppTarget;

    suppSize=size(supp,1);
    suppTargetSize=size(suppTarget,1);

    assembly=obj.assembly;
    U=assembly.node.current_U_mat;

    A=size(U);
    NodeNum=A(1);
    Uhis=zeros(increStep,NodeNum,3);

    fprintf('Loading Analysis Start');
    
    % Find the external forces that is currently applied on the structure
    currentAppliedForce=zeros(3*NodeNum,1);    
    for i=1:NodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = ...
            assembly.node.current_ext_force_mat(i,:);
    end  

    % Track the equilibrium
    for i=1:increStep
        step=1;
        R=1;
        lambda=i;
        fprintf('Icrement = %d\n',i);        

        
        while and(step<iterMax,R>tol)

            Usupp=zeros(3*NodeNum,1);
            for j=1:suppTargetSize
                TempNodeNum=suppTarget(j,1);
                Usupp(TempNodeNum*3-2)=U(TempNodeNum,1)-suppTarget(j,2)*i/increStep;
                Usupp(TempNodeNum*3-1)=U(TempNodeNum,2)-suppTarget(j,3)*i/increStep;
                Usupp(TempNodeNum*3)=U(TempNodeNum,3)-suppTarget(j,4)*i/increStep;
            end

            % find the internal force and stiffness of system
            [T,K]=assembly.Solve_FK(U);

            % calculate the unbalanced force
            unbalance=currentAppliedForce-T-K*Usupp; 
            
            [K,unbalance]=Mod_K_For_Supp(K,supp,unbalance);
            K=sparse(K);
                        
            dUtemp=(K\unbalance);
            for j=1:NodeNum
                U((j),:)=U((j),:)+dUtemp(3*j-2:3*j)';
            end
            for j=1:suppTargetSize
                TempNodeNum=suppTarget(j,1);
                U(TempNodeNum,1)=suppTarget(j,2)*i/increStep;
                U(TempNodeNum,2)=suppTarget(j,3)*i/increStep;
                U(TempNodeNum,3)=suppTarget(j,4)*i/increStep;
            end

            R=norm(dUtemp);
            fprintf('    Iteration = %d, R = %e\n',step,R);
            step=step+1;        
        end

        Uhis(i,:,:)=U;
    end  

    assembly.node.current_U_mat=U;
end