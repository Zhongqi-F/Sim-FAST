%% Plot the configuration of the model

function Plot_Shape_SprFlatNumber(obj)

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

barNum=size(assembly.bar.barConnect_Mat);
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
A=size(assembly.node.coordinates_Mat);
N=A(1);

sprIJK=obj.assembly.rotSprFlat.rotSprIJK_Mat;
sprNum=size(sprIJK);
sprNum=sprNum(1);


for i=1:sprNum
    x=node0(sprIJK(i,2),1);
    y=node0(sprIJK(i,2),2);
    z=node0(sprIJK(i,2),3);
    text(x,y,z,num2str(i),'Color','blue');
end

