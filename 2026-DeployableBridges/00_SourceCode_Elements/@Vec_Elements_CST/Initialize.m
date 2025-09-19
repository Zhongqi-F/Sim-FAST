function Initialize(obj,node)

    % There are multiple items that needs to be done to initialize
    % 1. We need to make sure that alpha2 and alpha 3 are not 90 degree.
    %    This is done by setting the largest angle to be alpha 1. 
    % 2. Calculate the L_mat.
    % 3. Solve the A_vec.


    numCST=size(obj.node_ijk_mat);
    numCST=numCST(1);

    obj.A_vec=zeros(numCST,1);

    for i=1:numCST
        % Identify the nodal coordinates
        n1=obj.node_ijk_mat(i,1);
        n2=obj.node_ijk_mat(i,2);
        n3=obj.node_ijk_mat(i,3);

        nodeIndex=[n1,n2,n3];

        % The three nodal coordinates of nodes
        X1=node.coordinates_mat(n1,:);
        X2=node.coordinates_mat(n2,:);
        X3=node.coordinates_mat(n3,:);

        % Find the each sector angle
        vtemp1=(X2-X1)/norm(X2-X1);
        vtemp2=(X3-X1)/norm(X3-X1);
        beta1=acos(dot(vtemp1,vtemp2));

        vtemp1=(X1-X2)/norm(X1-X2);
        vtemp2=(X3-X2)/norm(X3-X2);
        beta2=acos(dot(vtemp1,vtemp2));

        vtemp1=(X1-X3)/norm(X1-X3);
        vtemp2=(X2-X3)/norm(X2-X3);
        beta3=acos(dot(vtemp1,vtemp2));

        % Rank the sector angle by sizes
        betaVec=[beta1,beta2,beta3];
        [betaVec,index]=sort(betaVec);

        % Reorganize the node sequence
        obj.node_ijk_mat(i,1)=nodeIndex(index(3));
        obj.node_ijk_mat(i,2)=nodeIndex(index(2));
        obj.node_ijk_mat(i,3)=nodeIndex(index(1));
    end

    % Now that everything is organized in the right sequence, we can
    % calcualte the values we need for the elements. It is okay to not
    % optimize this part of the code, as we only use it once
    for i=1:numCST
        % Identify the nodal coordinates
        n1=obj.node_ijk_mat(i,1);
        n2=obj.node_ijk_mat(i,2);
        n3=obj.node_ijk_mat(i,3);

        % The three nodal coordinates of nodes
        X1=node.coordinates_mat(n1,:);
        X2=node.coordinates_mat(n2,:);
        X3=node.coordinates_mat(n3,:);

        L1=norm(X2-X3);
        L2=norm(X1-X3);
        L3=norm(X1-X2);

        % Store the Length into the storage matrix
        obj.L_mat(i,1)=L1;
        obj.L_mat(i,2)=L2;
        obj.L_mat(i,3)=L3;

        % Find the area of the triangle
        v1=X3-X1;
        v2=X2-X1;

        area=cross(v1,v2);
        area=norm(area)/2;

        % Store the are info in the storage vector
        obj.A_vec(i)=area;

    end
end

