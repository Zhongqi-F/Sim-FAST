
function [Tw,Kw]=SolveFK(obj,node,U)
    
        [Tw]=Wedge_GlobalForce(obj,node,U);
        [Kw]=Wedge_GlobalStiff(obj,node,U);

        % Here we solve the potential energy vector
        nodalCoordinates=node.coordinates_Mat;
        wedgeConnect=obj.wedgeConnect_Mat;
        L0ref=obj.L0ref_Mat;
        Kwedge=obj.Kwedge_Vec;
    
        wedgeNum=size(Kwedge);
        wedgeNum=wedgeNum(1);

        obj.currentStrainEnergy_Vec=zeros(wedgeNum,1);

        for i=1:wedgeNum

            node1=wedgeConnect(i,1);
            node2=wedgeConnect(i,2);
            node3=wedgeConnect(i,3);
            node4=wedgeConnect(i,4);
            node5=wedgeConnect(i,5);
            node6=wedgeConnect(i,6);
            
            X1=nodalCoordinates(node1,:)+U(node1,:);
            X2=nodalCoordinates(node2,:)+U(node2,:);
            X3=nodalCoordinates(node3,:)+U(node3,:);
            X4=nodalCoordinates(node4,:)+U(node4,:);
            X5=nodalCoordinates(node5,:)+U(node5,:);
            X6=nodalCoordinates(node6,:)+U(node6,:);

            X=[X1;X2;X3;X4;X5;X6];

            obj.currentStrainEnergy_Vec(1)=...
                obj.Wedge_Potential(X,Kwedge(i),L0ref(i,:));
        end

end