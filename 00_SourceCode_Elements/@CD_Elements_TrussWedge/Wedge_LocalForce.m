
function Flocal=Wedge_LocalForce(obj,X,E,A0ref,L0ref)

    % The wedge elemements has 6 nodes, the gradient should be a 18 by 1
    % vector. Local force is the numerical gradient of the bar potential.
    Flocal=zeros(18,1);
    
    % We use this small step delta to calculate numerical gradient
    delta=obj.delta;

    % for each element in the matrix, we do the following calculation
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
    
        % This is the central deference equation
        Flocal(i)=1/2/delta*(obj.Wedge_Potential(TempXfor,E,A0ref,L0ref)-...
            obj.Wedge_Potential(TempXback,E,A0ref,L0ref));

    end
end