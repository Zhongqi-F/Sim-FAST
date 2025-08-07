function [Tcontact]=Solve_Global_Force(obj,node,U)

    % retrieve information about the triangle
    tri_ijk=obj.tri_ijk_mat;

    % information about the node
    nodeNum=size(node.coordinates_mat);
    nodeNum=nodeNum(1);
    currentNode=node.coordinates_mat+U;

    delta=obj.delta;

    % information about contact group
    groupNum=obj.group_number;
    maxGroup=max(groupNum);
    groupInfo={};
    for i=1:maxGroup
        groupInfo{i}=find(groupNum==i);
    end

    % initialize
    Tcontact=zeros(nodeNum*3,1);

    % calculate contact forces
    for k=1:maxGroup
        for l=k+1:maxGroup
            for i=groupInfo{k}'
                for j=groupInfo{l}'
                    
                    nodeSet1=tri_ijk(i,:);
                    nodeSet2=tri_ijk(j,:);

                    xSet1=currentNode(nodeSet1,:);
                    xSet2=currentNode(nodeSet2,:);

                    % Find the distance
                    dist=obj.Solve_Distance(xSet1,xSet2);

                    if dist<obj.d0  
                        localF=obj.Solve_Local_Force(xSet1,xSet2,...
                            obj.d0,obj.k_contact,delta);
                    
                        % Add the force value to the appropriate 
                        % location in the global force vector
                        for t=1:3
                            nodeIndex=nodeSet1(t);
                            Tcontact((3*(nodeIndex-1)+1):(3*nodeIndex))=...
                                 Tcontact((3*(nodeIndex-1)+1):(3*nodeIndex))+...
                                 localF((3*(t-1)+1):(3*t));
                        end
                        for t=1:3
                            nodeIndex=nodeSet2(t);
                            Tcontact((3*(nodeIndex-1)+1):(3*nodeIndex))=...
                                 Tcontact((3*(nodeIndex-1)+1):(3*nodeIndex))+...
                                 localF((3*3+3*(t-1)+1):(3*3+3*t));
                        end
                    end

                end
            end
        end
    end

   
end