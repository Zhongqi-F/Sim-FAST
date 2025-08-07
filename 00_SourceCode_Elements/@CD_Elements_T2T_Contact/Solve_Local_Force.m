function Flocal=Solve_Local_Force(obj,xSet1,xSet2,d0,k,delta)

    % size of the local force vector is 3 times the total number of nodes.
    Flocal=zeros(3*6,1);
    
    % % We use this small step delta to calculate numerical gradient
    % delta=obj.delta;

    for i=1:3
        for j=1:3
            % forward step
            TempX1for=xSet1;
            TempX1for(i,j)=TempX1for(i,j)+delta;
            
            % backward step
            TempX1back=xSet1;
            TempX1back(i,j)=TempX1back(i,j)-delta;

            dfor=obj.Solve_Distance(TempX1for,xSet2);
            dback=obj.Solve_Distance(TempX1back,xSet2);

            % This is the central deference equation
            Flocal((i-1)*3+j)=1/2/delta*(obj.Potential(dfor,k,d0)-...
                obj.Potential(dback,k,d0));

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

            dfor=obj.Solve_Distance(xSet1,TempX2for);
            dback=obj.Solve_Distance(xSet1,TempX2back);

            % This is the central deference equation
            Flocal(3*3+(i-1)*3+j)=1/2/delta*(obj.Potential(dfor,k,d0)-...
                obj.Potential(dback,k,d0));    
           
        end
    end
end