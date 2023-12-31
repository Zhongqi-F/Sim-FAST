%% Plot the configuration of the model

function Plot_Shape_SprNumber(obj)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;

% Plot Dots
figure
view(View1,View2); 
set(gca,'DataAspectRatio',[1 1 1])
set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])

% The software support two ways to set up the plotting range
A=size(Vsize);
if A(1)==1    
    axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
else
    axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
end

wedgeConnect=assembly.wedge.wedgeConnect_Mat;
wedgeNum=size(wedgeConnect);
wedgeNum=wedgeNum(1);
node=assembly.node.coordinates_Mat;

for j=1:wedgeNum
    node1=node(wedgeConnect(j,1),:);
    node2=node(wedgeConnect(j,2),:);
    node3=node(wedgeConnect(j,3),:);
    node4=node(wedgeConnect(j,4),:);
    node5=node(wedgeConnect(j,5),:);
    node6=node(wedgeConnect(j,6),:);

    f=[1,2,3];
    v=[node1;node2;node3;];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3];
    v=[node4;node5;node6;];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node1;node2;node5;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node2;node3;node6;node5];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

    f=[1,2,3,4];
    v=[node1;node3;node6;node4];
    patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',.1)

end

barNum=size(assembly.bar.A_Vec);
barNum=barNum(1);
barConnect=assembly.bar.barConnect_Mat;

for j=1:barNum
    node1=assembly.node.coordinates_Mat(barConnect(j,1),:);
    node2=assembly.node.coordinates_Mat(barConnect(j,2),:);
    line([node1(1),node2(1)],...
         [node1(2),node2(2)],...
         [node1(3),node2(3)],'Color','k');
end

% Number Dots
node0=assembly.node.coordinates_Mat;

sprIJKL=obj.assembly.spr.sprIJKL_Mat;
sprNum=size(sprIJKL);
sprNum=sprNum(1);

for i=1:sprNum
    x=0.5*(node0(sprIJKL(i,2),1)+...
        node0(sprIJKL(i,3),1));
    y=0.5*(node0(sprIJKL(i,2),2)+...
        node0(sprIJKL(i,3),2));
    z=0.5*(node0(sprIJKL(i,2),3)+...
        node0(sprIJKL(i,3),3));
    text(x,y,z,num2str(i),'Color','blue');
end

