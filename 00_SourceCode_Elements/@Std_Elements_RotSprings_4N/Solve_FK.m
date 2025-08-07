function [Trs,Krs]=Solve_FK(obj,node,U)
    
        [Trs]=Solve_Global_Force(obj,node,U);
        [Krs]=Solve_Global_Stiff(obj,node,U);

        rotSprNum=size(obj.rot_spr_K_vec);
        rotSprNum=rotSprNum(1);

        for i=1:rotSprNum

            node1=obj.node_ijkl_mat(i,1);
            node2=obj.node_ijkl_mat(i,2);
            node3=obj.node_ijkl_mat(i,3);
            node4=obj.node_ijkl_mat(i,4);
            
            % The nodal cooridnates of the node after adding the deformation
            X1=node.coordinates_mat(node1,:)+U(node1,:);
            X2=node.coordinates_mat(node2,:)+U(node2,:);
            X3=node.coordinates_mat(node3,:)+U(node3,:);
            X4=node.coordinates_mat(node4,:)+U(node4,:);
            
            X=[X1;X2;X3;X4;];

            obj.theta_current_vec(i)=obj.Solve_Theta(X);
            obj.energy_current_vec(i)=obj.Potential(obj.theta_current_vec(i), ...
                obj.theta_stress_free_vec(i),obj.rot_spr_K_vec(i));
        end
        
end