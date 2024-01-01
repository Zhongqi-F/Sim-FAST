%% Initialize the solver
clear all;
clc;
close all;


%% Define the Geometry of origami
% Here we generate the geometry of the origami

% Define the nodal coordinate before meshing
R=50*10^(-3);
H=50*10^(-3);
theta=30/180*pi;
N=6;
M=2;

% Stiffness parameters of the structure
sprStiff=0.000001;
stiffFactor=1000;
barA=1*10^(-3)*10*10^(-3); 
barE=2*10^9; % Young's modulus

bar=Elements_Bars;
rotSpr=Elements_RotSprings;
node=Elements_Nodes;


%% Geometry of Kresling origami

alpha=2*pi/N;
for i=1:M+1

    for j=1:N
        node.coordinates_Mat=[node.coordinates_Mat;
            R*cos(j*alpha+i*theta),R*sin(j*alpha+i*theta),(i-1)*H];
    end
end

% Add the bar element that is horizontal
for i=1:M+1
    for j=1:N
        if j ~=N
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i-1)*N+j+1];
        else
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i-1)*N+1];
        end
    end
end

% Add those diagonal bars
for i=1:M
    for j=1:N
        if j ~=N
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i)*N+j];
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i)*N+j+1];
        else
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i)*N+j];
            bar.barConnect_Mat=[bar.barConnect_Mat;
                (i-1)*N+j,(i)*N+1];
        end
    end
end

% Initialize other parameters
barNum=size(bar.barConnect_Mat);
barNum=barNum(1);
bar.A_Vec=barA*ones(barNum,1);
bar.E_Vec=barE*ones(barNum,1);

%% Set up the rotational springs
% Diagonal rotational springs
for i=1:M
    for j=1:N
        if j ==1
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i-1)*N+N,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        elseif j~=N
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+j+1];
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i)*N+j,(i-1)*N+j,(i)*N+j+1,(i-1)*N+j+1];
        else
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i-1)*N+j-1,(i-1)*N+j,(i)*N+j,(i)*N+1];
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i)*N+j,(i-1)*N+j,(i)*N+1,(i-1)*N+1];
        end
    end
end

for i=1:M-1
    for j=1:N
        if j~=N
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i-1)*N+j,i*N+j,i*N+j+1,(i+1)*N+j+1];
        else
            rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;
                (i-1)*N+j,i*N+j,i*N+1,(i+1)*N+1];
        end
    end
end

rotSpr.rotSprK_Vec=sprStiff*ones(30,1);



%% Initialize assembly
assembly=Assembly_Origami();
assembly.node=node;
assembly.bar=bar;
assembly.rotSpr=rotSpr;

assembly.InitializeAssembly()

%% Plot for investigation
plots=Plot_Origami();
plots.displayRange=0.3;
plots.displayRangeRatio=1;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()
plots.Plot_Shape_SprNumber()



%% Setup the loading controller
dc=Solver_DC;
dc.assembly=assembly;
dc.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;
         4,1,1,1;
         5,1,1,1;
         6,1,1,1];

force=10;

dc.selectedRefDisp=[13,3];

dc.load=[13,0,0,-force;
         14,0,0,-force;
         15,0,0,-force;
         16,0,0,-force;
         17,0,0,-force;
         18,0,0,-force;];

dc.increStep=50;
dc.tol=10^-5;
dc.iterMax=50;

Uhis=dc.Solve();

plots.displayRange=0.15;
plots.fileName='Kresling.gif';
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)));
plots.Plot_DeformedHis(Uhis);

%% Find the reaction force and loading results

forceHis=zeros(dc.increStep,1);
UrefHis=zeros(dc.increStep,1);

for i=1:dc.increStep
    [F,K]=assembly.SolveFK(squeeze(Uhis(i,:,:)));
    UrefHis(i)=Uhis(i,dc.selectedRefDisp(1),dc.selectedRefDisp(2));
    forceHis(i)=F(dc.selectedRefDisp(2)+(dc.selectedRefDisp(1))*3);
end

figure
plot(UrefHis,forceHis)
xlabel('Z deformation of top node (m)') 
ylabel('Applied Force (N)') 


EnergyHis=zeros(dc.increStep,2);

for i=1:dc.increStep
    % rotational spring element
    rotSpr.CalcStrainEnergy(node,squeeze(Uhis(i,:,:)));
    EnergyHis(i,1)=sum(rotSpr.strainEnergy_Vec);

    % bar element
    bar.CalcStrainEnergy(node,squeeze(Uhis(i,:,:)));
    EnergyHis(i,2)=sum(bar.currentStrainEnergy_Vec);

end

figure
hold on
plot(UrefHis,EnergyHis(:,1))
plot(UrefHis,EnergyHis(:,2))
xlabel('Z deformation of top node (m)') 
ylabel('Strain Energy (J)')
legend('rotational springs','bars')