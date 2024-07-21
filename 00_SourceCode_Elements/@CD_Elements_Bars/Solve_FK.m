function [Tbar,Kbar]=Solve_FK(obj,node,U)
    
        [Tbar]=Solve_Global_Force(obj,node,U);
        [Kbar]=Solve_Global_Stiff(obj,node,U);

        % Here we solve the potential energy vector
        % Retrieve information about the bars
        nodalCoordinates=node.coordinates_mat;
        barConnect=obj.node_ij_mat;
        L0=obj.L0_vec;
        barA=obj.A_vec;
        barE=obj.E_vec;
    
        barNum=size(barA);
        barNum=barNum(1);

        obj.energy_current_vec=zeros(barNum,1);
        obj.strain_current_vec=zeros(barNum,1);

        % Compute the potential energy of each bar
        % This is used for output data.
        for i=1:barNum

            node1=barConnect(i,1);
            node2=barConnect(i,2);
            
            X1=nodalCoordinates(node1,:)+U(node1,:);
            X2=nodalCoordinates(node2,:)+U(node2,:);

            obj.strain_current_vec(i)=...
                (norm(X1-X2)-L0(i))/L0(i);

            obj.energy_current_vec(i)=...
                1/2*barE(i)*barA(i)/L0(i)*(norm(X1-X2)-L0(i))^2;
        end

end