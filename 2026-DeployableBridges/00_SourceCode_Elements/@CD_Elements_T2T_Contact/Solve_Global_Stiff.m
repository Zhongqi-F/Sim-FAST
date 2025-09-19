function [Kcontact]=Solve_Global_Stiff(obj,node,U)
    
    % retrieve information about the contacting pair
    tri_ijk=obj.tri_ijk_mat;

    % information about the node
    nodeNum=size(node.coordinates_mat);
    nodeNum=nodeNum(1);
    currentNode=node.coordinates_mat+U;

    % information about contact group
    groupNum=obj.group_number;
    maxGroup=max(groupNum);
    groupInfo={};
    for i=1:maxGroup
        groupInfo{i}=find(groupNum==i);
    end

    % initialize
    Kcontact=zeros(nodeNum*3,nodeNum*3);

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

                    if dist < obj.d0
                        localK=obj.Solve_Local_Stiff(xSet1,xSet2,obj.d0,obj.k_contact);
            
                        % Add the stiffness value to the appropriate location in the 
                        % global stiffness matrix. The local stiffness can be devided 
                        % into four different blocks, the following code addresses the 
                        % four blocks separately.
            
                        % block 1 (objec 1 : object 1)
                        for t=1:3
                            for p=1:3
                                nodeIndex1=nodeSet1(t);
                                nodeIndex2=nodeSet1(p);
                                Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))=...
                                    Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))+...
                                    localK((3*(t-1)+1):(3*t),(3*(p-1)+1):(3*p));
                            end
                        end
                        
                        % block 1 (objec 1 : object 2)
                        for t=1:3
                            for p=1:3
                                nodeIndex1=nodeSet1(t);
                                nodeIndex2=nodeSet2(p);
                                Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))=...
                                    Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))+...
                                    localK((3*(t-1)+1):(3*t),...
                                    (3*3+3*(p-1)+1):(3*3+3*p));
                            end
                        end
            
                        % block 1 (objec 2 : object 1)
                        for t=1:3
                            for p=1:3
                                nodeIndex1=nodeSet2(t);
                                nodeIndex2=nodeSet1(p);
                                Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))=...
                                    Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))+...
                                    localK((3*3+3*(t-1)+1):(3*3+3*t),...
                                    (3*(p-1)+1):(3*p));
                            end
                        end
                        
                        % block 1 (objec 2 : object 2)
                        for t=1:3
                            for p=1:3
                                nodeIndex1=nodeSet2(t);
                                nodeIndex2=nodeSet2(p);
                                Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))=...
                                    Kcontact((3*(nodeIndex1-1)+1):(3*nodeIndex1), ...
                                    (3*(nodeIndex2-1)+1):(3*nodeIndex2))+...
                                    localK((3*3+3*(t-1)+1):(3*3+3*t),...
                                    (3*3+3*(p-1)+1):(3*3+3*p));
                            end
                        end

                    end
                end
            end
        end
    end    
end