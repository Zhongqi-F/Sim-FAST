function Klocal=Solve_Local_Stiff(obj,xSet1,xSet2,d0,k)

    % The dimension of the stiffness matrix is 3 times the total number of
    % nodes within the two body
    Klocal=zeros(6*3);
        
    % We use this small step delta to calculate numerical Hessian
    delta=obj.delta*10;

    for i=1:3
        for j=1:3
            % forward step, where the coodinates is moved in x direction 
            % by delta value for node 1
            TempX1for=xSet1;
            TempX1for(i,j)=TempX1for(i,j)+delta;
            
            % backward step, where the coodinates is moved in x direction
            % by -delta value for node 1
            TempX1back=xSet1;
            TempX1back(i,j)=TempX1back(i,j)-delta;
    
            % This is the central difference step
            Klocal((i-1)*3+j,:)=1/2/delta*(...
                obj.Solve_Local_Force(TempX1for,xSet2,d0,k,delta/10)...
                -obj.Solve_Local_Force(TempX1back,xSet2,d0,k,delta/10));

        end
    end

    for i=1:3
        for j=1:3
            % forward step
            TempX2for=xSet2;
            TempX2for(i,j)=TempX2for(i,j)+delta;
            
            % backward step
            TempX2back=xSet2;
            TempX2back(i,j)=TempX2back(i,j)-delta;
    
            % This is the central difference step
            Klocal(3*3+(i-1)*3+j,:)=1/2/delta*(...
                obj.Solve_Local_Force(xSet1,TempX2for,d0,k,delta/10)...
                -obj.Solve_Local_Force(xSet1,TempX2back,d0,k,delta/10));

        end
    end
end