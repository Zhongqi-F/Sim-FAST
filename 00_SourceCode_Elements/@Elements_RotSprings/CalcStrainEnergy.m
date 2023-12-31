%% Solve the strain energy and current fold angle of rotational springs
% This function calcualtes the strain energy and current folding angle 
% of the rotational springs. This function is used after solving the 
% loading or folding process. We use this function to compute responses
% (strain energy) that are of interest to us. 

function CalcStrainEnergy(obj,node,U)
    
        % find the number of elements
        rotSprNum=size(obj.rotSprK_Vec);
        rotSprNum=rotSprNum(1);

        for i=1:rotSprNum

            node1=obj.rotSprIJKL_Mat(i,1);
            node2=obj.rotSprIJKL_Mat(i,2);
            node3=obj.rotSprIJKL_Mat(i,3);
            node4=obj.rotSprIJKL_Mat(i,4);
            
            % The nodal cooridnates of the node after adding the deformation
            X1=node.coordinates_Mat(node1,:)+U(node1,:);
            X2=node.coordinates_Mat(node2,:)+U(node2,:);
            X3=node.coordinates_Mat(node3,:)+U(node3,:);
            X4=node.coordinates_Mat(node4,:)+U(node4,:);
            
            X=[X1;X2;X3;X4;];

            obj.theta_Current_Vec(i)=obj.Theta(X);
            obj.strainEnergy_Vec(i)=obj.Potential(obj.theta_Current_Vec(i), ...
                obj.theta_StressFree_Vec(i),obj.rotSprK_Vec(i));
        end

end