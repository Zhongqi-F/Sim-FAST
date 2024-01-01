%% This function initialize the assembly system

function AddCompliantCrease(obj,LeftNode,RightNode,LeftAnchor,RightAnchor,t,W,E,G)

    %% Add the three additional nodes
    coordinates=obj.node.coordinates_Mat;

    newX1=(coordinates(LeftNode(1),:)+...
        coordinates(RightNode(1),:))/2;
    newX2=(coordinates(LeftNode(2),:)+...
        coordinates(RightNode(2),:))/2;
    newX3=(newX1+newX2)/2;

    nodeNum=size(coordinates);
    nodeNum=nodeNum(1);
    newN1=nodeNum+1;
    newN2=nodeNum+2;
    newN3=nodeNum+3;
    
    coordinates=[coordinates;newX1];
    coordinates=[coordinates;newX2];
    coordinates=[coordinates;newX3];

    obj.node.coordinates_Mat=coordinates;

    L=norm(newX1-newX3);

    %% Add the bars
    % four longitudinal bar
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;LeftNode(1),newN1];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;RightNode(1),newN1];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;LeftNode(2),newN2];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;RightNode(2),newN2];

    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];

    % two horizontal bar
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;newN1,newN3];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;newN2,newN3];

    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t/2];

    % four diagonal bar
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;LeftNode(1),newN3];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;RightNode(1),newN3];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;LeftNode(2),newN3];
    obj.bar.barConnect_Mat=[obj.bar.barConnect_Mat;RightNode(2),newN3];

    obj.bar.A_Vec=[obj.bar.A_Vec;L*t*G/E/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t*G/E/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t*G/E/2];
    obj.bar.A_Vec=[obj.bar.A_Vec;L*t*G/E/2];

    obj.bar.E_Vec=[obj.bar.E_Vec;E*ones(10,1)];

    %% Add the rotational springs

    sprStiff=3*E*L*(t^3)/W/12;

    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;LeftAnchor,LeftNode(1),LeftNode(2),newN3];
    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;LeftNode(2),LeftNode(1),newN3,newN1];
    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;LeftNode(1),LeftNode(2),newN3,newN2];

    obj.rotSpr.rotSprK_Vec=[obj.rotSpr.rotSprK_Vec;sprStiff*ones(3,1)];

    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;LeftNode(1),newN3,newN1,RightNode(1)];
    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;LeftNode(2),newN3,newN2,RightNode(2)];

    obj.rotSpr.rotSprK_Vec=[obj.rotSpr.rotSprK_Vec;4*sprStiff*ones(2,1)];

    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;RightNode(2),RightNode(1),newN3,newN1];
    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;RightNode(1),RightNode(2),newN3,newN2];
    obj.rotSpr.rotSprIJKL_Mat=[obj.rotSpr.rotSprIJKL_Mat;RightAnchor,RightNode(1),RightNode(2),newN3];

    obj.rotSpr.rotSprK_Vec=[obj.rotSpr.rotSprK_Vec;sprStiff*ones(3,1)];
    
    obj.rotSpr.theta_StressFree_Vec=[obj.rotSpr.theta_StressFree_Vec;pi*ones(8,1)];

end