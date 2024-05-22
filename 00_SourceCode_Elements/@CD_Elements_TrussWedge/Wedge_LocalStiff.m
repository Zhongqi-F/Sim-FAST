
function Klocal=Wedge_LocalStiff(obj,X,E,A0ref,L0ref)

    % The wedge elemements has 6 nodes, the hessian should be a 18 by 18
    % matrix. Local stiffness is the numerical hessian of the potential.
    Klocal=zeros(18,18);
    
    % We use this small step delta to calculate numerical Hessian
    delta=obj.delta;

    
    for i=1:18
        % find the index
        if mod(i,3)==0
            index1=3;
            index2=i/3;
        else
            index1=mod(i,3);
            index2=floor(i/3)+1;
        end

        TempXfor=X;
        TempXfor(index2,index1)=TempXfor(index2,index1)+delta;
        
        % backward step, where the coodinates is moved in x direction by -delta 
        % value for node 1
        TempXback=X;
        TempXback(index2,index1)=TempXback(index2,index1)-delta;
    
        % This is the central deference equation. Please note that we can make 
        % use of the force function for this calculation. 
        Klocal(i,:)=1/2/delta*(obj.Wedge_LocalForce(TempXfor,E,A0ref,L0ref)-...
            obj.Wedge_LocalForce(TempXback,E,A0ref,L0ref));


    end    
end