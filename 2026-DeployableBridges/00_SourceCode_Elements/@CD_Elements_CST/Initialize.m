function Initialize(obj,node)

    numCST=size(obj.node_ijk_mat);
    numCST=numCST(1);

    for i=1:numCST

        n1=obj.node_ijk_mat(i,1);
        n2=obj.node_ijk_mat(i,2);
        n3=obj.node_ijk_mat(i,3);

        nodeIndex=[n1,n2,n3];

        x1=node.coordinates_mat(n1,:);
        x2=node.coordinates_mat(n2,:);
        x3=node.coordinates_mat(n3,:);

        vtemp1=(x2-x1)/norm(x2-x1);
        vtemp2=(x3-x1)/norm(x3-x1);
        alpha1=acos(dot(vtemp1,vtemp2));

        vtemp1=(x1-x2)/norm(x1-x2);
        vtemp2=(x3-x2)/norm(x3-x2);
        alpha2=acos(dot(vtemp1,vtemp2));

        vtemp1=(x1-x3)/norm(x1-x3);
        vtemp2=(x2-x3)/norm(x2-x3);
        alpha3=acos(dot(vtemp1,vtemp2));

        alphaVec=[alpha1,alpha2,alpha3];
        [alphaVec,index]=sort(alphaVec);

        obj.node_ijk_mat(i,1)=nodeIndex(index(3));
        obj.node_ijk_mat(i,2)=nodeIndex(index(2));
        obj.node_ijk_mat(i,3)=nodeIndex(index(1));

    end
end

