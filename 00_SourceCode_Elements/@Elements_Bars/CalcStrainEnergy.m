%% Solve the stored strain energy

function CalcStrainEnergy(obj,node,U)

        % Here we solve the potential energy vector
        % Retrieve information about the bars
        nodalCoordinates=node.coordinates_Mat;
        barConnect=obj.barConnect_Mat;
        L0=obj.L0_Vec;
        barA=obj.A_Vec;
        barE=obj.E_Vec;
    
        barNum=size(barA);
        barNum=barNum(1);

        obj.currentStrainEnergy_Vec=zeros(barNum,1);

        % Compute the potential energy of each bar
        % This is used for output data.
        for i=1:barNum

            node1=barConnect(i,1);
            node2=barConnect(i,2);
            
            X1=nodalCoordinates(node1,:)+U(node1,:);
            X2=nodalCoordinates(node2,:)+U(node2,:);

            obj.currentStrainEnergy_Vec(1)=...
                1/2*barE(i)*barA(i)/L0(i)*(norm(X1-X2)-L0(i))^2;
        end

end