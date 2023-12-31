%% Nonlinear solver for folding (Assembling)

function [Uhis]=Solve(obj)
    
    increStep=obj.increStep;
    tol=obj.tol;
    iterMax=obj.iterMax;    
    supp=obj.supp;  

    assembly=obj.assembly;
    U=assembly.node.currentU_Mat;

    A=size(U);
    NodeNum=A(1);
    Uhis=zeros(increStep,NodeNum,3);

    % Find the external forces that is currently applied on the structure
    currentAppliedForce=zeros(3*NodeNum,1);    
    for i=1:NodeNum
        currentAppliedForce(3*(i-1)+1:3*i) = assembly.node.currentExtForce_Mat(i,:);
    end  


   
    % Analysis start
    fprintf('Self Assemble Analysis Start');
    count=1;

    % find the zero strain angle, before and after the analysis
    L0_before=assembly.bar.L0_Vec;
    L0_after=obj.targetL0;
    
        
    for i=1:increStep
        
        L0_current=i/increStep*L0_after+...
            (1-i/increStep)*L0_before;

        fprintf('Icrement = %d\n',count);
        step=1;
        R=1;
             
        while and(step<iterMax,R>tol)
            
            assembly.bar.L0_Vec=...
                L0_current;
            [T,K]=assembly.SolveFK(U);

            % calculate the unbalanced force
            unbalance=currentAppliedForce-T; 
            
            [K,unbalance]=obj.ModKforSupp(K,supp,unbalance);
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
end

