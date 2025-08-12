function Add_Triangle_Panel(obj,n1,n2,n3,n4,n5,n6,E,t,v)

        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n1,n2,n3];
        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n4,n5,n6];

        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n1,n2,n5];
        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n1,n4,n5];

        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n2,n3,n6];
        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n2,n5,n6];

        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n3,n1,n4];
        obj.cst.node_ijk_mat=[obj.cst.node_ijk_mat;n3,n6,n4];

        obj.cst.E_vec=[obj.cst.E_vec;E*ones(8,1)];
        obj.cst.t_vec=[obj.cst.t_vec;t*ones(8,1)];
        obj.cst.v_vec=[obj.cst.v_vec;v*ones(8,1)];

end

