clear;clc;close all;

%% Set up the geometry of the knitting structure
node=Elements_Nodes();
bar=Elements_Bars();
rotSpr=Elements_RotSprings();
rotSprFlat=Elements_RotSprings_Flat();

% number of cells
M=3;
N=3;

% grid dimension
L=0.05;
offset=0.01;

% Stiffness parameters
barArea=0.00001;
barE=2*10^9;
bendStiff=1;
twistStiff=100;

% bar area in the z-direction
zBarFactor=0.0001;

% Applied load
force=1;

% The number of node after adding structures in y direction
NodeNumAfterY=(M)*(N+2);

% Add structure in the y direction
for i=1:M
    node.coordinates_Mat=[node.coordinates_Mat;
        0, (i-1)*L+L/2, 0]; 
    for j=1:N
        node.coordinates_Mat=[node.coordinates_Mat;
            (j-1)*L+L/2, (i-1)*L+L/2, ((-1)^j)*((-1)^i)*offset/2]; 
        bar.barConnect_Mat=[bar.barConnect_Mat;
            (i-1)*(N+2)+j,(i-1)*(N+2)+j+1];
        % Add rotational springs (3 node flat version)
        rotSprFlat.rotSprIJK_Mat=[rotSprFlat.rotSprIJK_Mat;
            (i-1)*(N+2)+j,(i-1)*(N+2)+j+1,NodeNumAfterY+(j-1)*(M+2)+i+1];
        rotSprFlat.rotSprIJK_Mat=[rotSprFlat.rotSprIJK_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i+1,(i-1)*(N+2)+j+1,(i-1)*(N+2)+j+2];
    end
    node.coordinates_Mat=[node.coordinates_Mat;
        N*L, (i-1)*L+L/2, 0]; 
    bar.barConnect_Mat=[bar.barConnect_Mat;
            i*(N+2)-1,i*(N+2)];
end

% Add structure in the x direction
for j=1:N
    node.coordinates_Mat=[node.coordinates_Mat;
        (j-1)*L+L/2, 0, 0]; 
    for i=1:M
        node.coordinates_Mat=[node.coordinates_Mat;
            (j-1)*L+L/2, (i-1)*L+L/2, -((-1)^j)*((-1)^i)*offset/2];
        % Add the bar element
        bar.barConnect_Mat=[bar.barConnect_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i,NodeNumAfterY+(j-1)*(M+2)+i+1];
        % Add rotational springs (3 node flat version)
        rotSprFlat.rotSprIJK_Mat=[rotSprFlat.rotSprIJK_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i,NodeNumAfterY+(j-1)*(M+2)+i+1,...
            (i-1)*(N+2)+j+1];
        rotSprFlat.rotSprIJK_Mat=[rotSprFlat.rotSprIJK_Mat;
            (i-1)*(N+2)+j+1,NodeNumAfterY+(j-1)*(M+2)+i+1,...
            NodeNumAfterY+(j-1)*(M+2)+i+2];
    end
    node.coordinates_Mat=[node.coordinates_Mat;
        (j-1)*L+L/2, M*L, 0];
    bar.barConnect_Mat=[bar.barConnect_Mat;
            NodeNumAfterY+j*(M+2)-1,NodeNumAfterY+j*(M+2)];
end

% Add bar in the z direction
for i=1:M
    for j=1:N
        bar.barConnect_Mat=[bar.barConnect_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i+1,(i-1)*(N+2)+j+1];
        
    end
end

% Add rotational springs (4 nodes)
for i=1:M
    for j=1:N
        rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
            (i-1)*(N+2)+j, NodeNumAfterY+(j-1)*(M+2)+i+1, ...
            (i-1)*(N+2)+j+1, NodeNumAfterY+(j-1)*(M+2)+i;];
        rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i, NodeNumAfterY+(j-1)*(M+2)+i+1, ...
            (i-1)*(N+2)+j+1, (i-1)*(N+2)+j+2];
        rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
            (i-1)*(N+2)+j+2, NodeNumAfterY+(j-1)*(M+2)+i+1, ...
            (i-1)*(N+2)+j+1, NodeNumAfterY+(j-1)*(M+2)+i+2];
        rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
            NodeNumAfterY+(j-1)*(M+2)+i+2, NodeNumAfterY+(j-1)*(M+2)+i+1, ...
            (i-1)*(N+2)+j+1, (i-1)*(N+2)+j];
    end
end


%% Set up the stiffness parameters
barNum=size(bar.barConnect_Mat);
barNum=barNum(1);
rotSprNum=size(rotSpr.rotSprIJKL_Mat);
rotSprNum=rotSprNum(1);
rotSprFlatNum=size(rotSprFlat.rotSprIJK_Mat);
rotSprFlatNum=rotSprFlatNum(1);


bar.A_Vec=barArea*ones(barNum,1);
bar.E_Vec=barE*ones(barNum,1);

verticalBarStart=(M+1)*N+(N+1)*M+1;
bar.A_Vec(verticalBarStart:end)=bar.A_Vec(verticalBarStart:end)/zBarFactor;

rotSpr.rotSprK_Vec=twistStiff*ones(rotSprNum,1);
rotSprFlat.rotSprK_Vec=bendStiff*ones(rotSprFlatNum,1);


%% Initialize assembly
assembly=Assembly_Knit();
assembly.node=node;
assembly.bar=bar;
assembly.rotSpr=rotSpr;
assembly.rotSprFlat=rotSprFlat;

assembly.InitializeAssembly()

%% Plot for investigation
plots=Plot_Knit();
plots.displayRange=0.15;
plots.displayRangeRatio=0.2;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()
plots.Plot_Shape_SprNumber()
plots.Plot_Shape_SprFlatNumber()


%% Check the mode
delta=10^-6;
rotSprFlat.delta=delta;
rotSpr.delta=delta;
bar.delta=delta;

nodeNum=(M+2)*N+(N+2)*M;
[F,K]=assembly.SolveFK(zeros(nodeNum,3));

[U,V]=eigs(K,15,'smallestabs');
D=diag(V)



%% Set up loading

nr=Solver_NR_Loading();
nr.assembly=assembly;


suppIndex=[1:(N+2):((N+2)*(M-1)+1),...
    (N+2):(N+2):(N+2)*M,];
suppIndex2=(NodeNumAfterY+1):(M+2):(NodeNumAfterY+(M+2)*(N-1)+1);

loadIndex=(NodeNumAfterY+M+2):(M+2):(NodeNumAfterY+(M+2)*N);

nr.supp=[suppIndex',ones((N)*2,1),zeros((N)*2,1),ones((N)*2,1)];
nr.supp=[nr.supp;
    suppIndex2',ones(M,3)];

nr.load=[loadIndex',zeros(M,1),force*ones(M,1),zeros(M,1)];

nr.increStep=40;
nr.tol=10^-5;

Uhis=nr.Solve();
Ufinal=squeeze(Uhis(nr.increStep,:,:));
plots.Plot_DeformedShape(Ufinal);

