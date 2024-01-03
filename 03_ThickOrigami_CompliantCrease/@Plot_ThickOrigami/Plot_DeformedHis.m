%% plot the deformation history of the simulated results.


function Plot_DeformedHis(obj,Uhis)

View1=obj.viewAngle1;
View2=obj.viewAngle2;
Vsize=obj.displayRange;
Vratio=obj.displayRangeRatio;

assembly=obj.assembly;
undeformedNode=assembly.node.coordinates_Mat;

pauseTime=obj.holdTime;
filename=obj.fileName;

h=figure;

set(gcf, 'color', 'white');
set(gcf,'position',[obj.x0,obj.y0,obj.width,obj.height])
  
B=size(Uhis);
Incre=B(1);

for i=1:Incre
    clf
    view(View1,View2); 
    set(gca,'DataAspectRatio',[1 1 1])
    
    A=size(Vsize);
    if A(1)==1    
        axis([-Vratio*Vsize Vsize -Vratio*Vsize Vsize -Vratio*Vsize Vsize])
    else
        axis([Vsize(1) Vsize(2) Vsize(3) Vsize(4) Vsize(5) Vsize(6)])
    end
    
    tempU=squeeze(Uhis(i,:,:));
    deformNode=undeformedNode+tempU;

    barConnect=assembly.bar.barConnect_Mat;
    barNum=size(barConnect);
    barNum=barNum(1);
    for j=1:barNum
        node1=deformNode(barConnect(j,1),:);
        node2=deformNode(barConnect(j,2),:);
        line([node1(1),node2(1)],...
             [node1(2),node2(2)],...
             [node1(3),node2(3)],'Color','k');
    end
    
    wedgeConnect=assembly.wedge.wedgeConnect_Mat;
    wedgeNum=size(wedgeConnect);
    wedgeNum=wedgeNum(1);

    for j=1:wedgeNum
        node1=deformNode(wedgeConnect(j,1),:);
        node2=deformNode(wedgeConnect(j,2),:);
        node3=deformNode(wedgeConnect(j,3),:);
        node4=deformNode(wedgeConnect(j,4),:);
        node5=deformNode(wedgeConnect(j,5),:);
        node6=deformNode(wedgeConnect(j,6),:);
    
        f=[1,2,3];
        v=[node1;node2;node3;];
        patch('Faces',f,'Vertices',v,'FaceColor','yellow')
    
        f=[1,2,3];
        v=[node4;node5;node6;];
        patch('Faces',f,'Vertices',v,'FaceColor','yellow')
    
        f=[1,2,3,4];
        v=[node1;node2;node5;node4];
        patch('Faces',f,'Vertices',v,'FaceColor','yellow')
    
        f=[1,2,3,4];
        v=[node2;node3;node6;node5];
        patch('Faces',f,'Vertices',v,'FaceColor','yellow')
    
        f=[1,2,3,4];
        v=[node1;node3;node6;node4];
        patch('Faces',f,'Vertices',v,'FaceColor','yellow')    
    
    end

    pause(pauseTime);  
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);

    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', pauseTime); 
    end 
end
