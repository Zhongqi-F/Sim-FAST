%% Initialize the solver
clear all;
clc;
close all;
tic


%% Define the Geometry of origami
% This section of code is used to generate the geometry of the origami
% pattern before meshing; 

% Definition of Geometry
L1=300*10^-6; % Panel side overhang
L2=100*10^-6; % Beam width
L3=600*10^-6; % Panel middle width
L4=400*10^-6; % Panel length
L5=440*10^-6; % Beam Length

DampingFactorA=100;
DampingFactorB=0;
ExciatationAngle=0.25;
totalSweepTime=0.8;

bar=Elements_Bars;
rotSpr=Elements_RotSprings;
node=Elements_Nodes;

%% Add node of the micro-mirror
% Here we define the nodal coordinates
node.coordinates_Mat(1,:)=[-L5/3,L3/2+L2,0];
node.coordinates_Mat(2,:)=[0,L3/2+L2,0];
node.coordinates_Mat(3,:)=[L5/3,L3/2+L2,0];
node.coordinates_Mat(4,:)=[L5*2/3,L3/2+L2,0];
node.coordinates_Mat(5,:)=[L5,L3/2+L2,0];

node.coordinates_Mat(6,:)=[0,L3/2,0];
node.coordinates_Mat(7,:)=[L5/3,L3/2,0];
node.coordinates_Mat(8,:)=[L5*2/3,L3/2,0];
node.coordinates_Mat(9,:)=[L5,L3/2,0];

node.coordinates_Mat(10,:)=[-L5/3,-L3/2-L2,0];
node.coordinates_Mat(11,:)=[0,-L3/2-L2,0];
node.coordinates_Mat(12,:)=[L5/3,-L3/2-L2,0];
node.coordinates_Mat(13,:)=[L5*2/3,-L3/2-L2,0];
node.coordinates_Mat(14,:)=[L5,-L3/2-L2,0];

node.coordinates_Mat(15,:)=[0,-L3/2,0];
node.coordinates_Mat(16,:)=[L5/3,-L3/2,0];
node.coordinates_Mat(17,:)=[L5*2/3,-L3/2,0];
node.coordinates_Mat(18,:)=[L5,-L3/2,0];

node.coordinates_Mat(19,:)=[L5,L3/2+L2+L1,0];
node.coordinates_Mat(20,:)=[L5,-L3/2-L2-L1,0];

node.coordinates_Mat(21,:)=[L5+L4/2,0,0];

node.coordinates_Mat(22,:)=[L5+L4,L3/2+L2+L1,0];
node.coordinates_Mat(23,:)=[L5+L4,L3/2,0];
node.coordinates_Mat(24,:)=[L5+L4,-L3/2,0];
node.coordinates_Mat(25,:)=[L5+L4,-L3/2-L2-L1,0];


%% solve for the area of bars 
% Layer information of the PZT beam
tPtTop=0.1*10^-6;
tPZT=1*10^-6;
tPtBTM=0.15*10^-6;
tSiO2=1*10^-6;

tAl=1*10^-6;

ESiO2=66*10^9;
EPt=168*10^9;
EPZT=70*10^9;
EAl=70*10^9;

Eeq=100*10^9;
% Assume an equivalent Young's modulus of 100 GPa

tBeamTotal=tSiO2+tPtTop+tPZT+tPtBTM;
EAbeam=L2*(tSiO2*ESiO2+tPtTop*EPt+tPZT*EPZT+tPtBTM*EPt);

% Assume material is homogenized with 100 GPa
AbeamLong=EAbeam/Eeq/2;
AbeamDiag=3*AbeamLong;
AbeamHor=3*AbeamLong;
Apanel=5*AbeamLong;


%% Solve for the rotational springs stiffness for bending 

AreaMomentBeam=(tSiO2*L2*(tSiO2/2)*ESiO2)+...
    (tPtBTM*L2*(tPtBTM/2+tSiO2)*EPt)+...
    (tPZT*L2*(tPZT/2+tPtBTM+tSiO2)*EPZT)+...
    (tPtTop*L2*(tPtTop/2+tPZT+tPtBTM+tSiO2)*EPt);
c=AreaMomentBeam/EAbeam;

EIcomp=L2*(ESiO2*(1/12*tSiO2^3+tSiO2*(c-tSiO2/2)^2)+...
    EPt*(1/12*tPtBTM^3+tPtBTM*(c-tPtBTM/2-tSiO2)^2)+...
    EPZT*(1/12*tPZT^3+tPZT*(c-tPZT/2-tPtBTM-tSiO2)^2)+...
    EPt*(1/12*tPtTop^3+tPtTop*(c-tPtTop/2-tPZT-tPtBTM-tSiO2)^2));


% tPary=1*10^-6;
% EPary=2*10^9;
% EAbeamPary=L2*(tSiO2*ESiO2+tPtTop*EPt+tPZT*EPZT+tPtBTM*EPt+tPary*EPary);
% AreaMomentBeamPary=(tSiO2*L2*(tSiO2/2)*ESiO2)+...
%     (tPtBTM*L2*(tPtBTM/2+tSiO2)*EPt)+...
%     (tPZT*L2*(tPZT/2+tPtBTM+tSiO2)*EPZT)+...
%     (tPtTop*L2*(tPtTop/2+tPZT+tPtBTM+tSiO2)*EPt)+...
%     (tPary*L2*(tPary/2+tPtTop+tPZT+tPtBTM+tSiO2)*EPary);
% c=AreaMomentBeamPary/EAbeamPary;
% EIcompPary=L2*(ESiO2*(1/12*tSiO2^3+tSiO2*(c-tSiO2/2)^2)+...
%     EPt*(1/12*tPtBTM^3+tPtBTM*(c-tPtBTM/2-tSiO2)^2)+...
%     EPZT*(1/12*tPZT^3+tPZT*(c-tPZT/2-tPtBTM-tSiO2)^2)+...
%     EPt*(1/12*tPtTop^3+tPtTop*(c-tPtTop/2-tPZT-tPtBTM-tSiO2)^2)+...
%     EPary*(1/12*tPary^3+tPary*(c-tPary/2-tPtTop-tPZT-tPtBTM-tSiO2)^2));

SprBeam=4*EIcomp/L5;
SprBeamDiag=3*SprBeam;
SprPanel=60*SprBeam;



%% Define longitudinal bar elements
bar.barConnect_Mat=[1,2];
bar.barConnect_Mat=[bar.barConnect_Mat; 2,3];
bar.barConnect_Mat=[bar.barConnect_Mat; 3,4];
bar.barConnect_Mat=[bar.barConnect_Mat; 4,5];

bar.barConnect_Mat=[bar.barConnect_Mat; 6,7];
bar.barConnect_Mat=[bar.barConnect_Mat; 7,8];
bar.barConnect_Mat=[bar.barConnect_Mat; 8,9];

bar.barConnect_Mat=[bar.barConnect_Mat; 10,11];
bar.barConnect_Mat=[bar.barConnect_Mat; 11,12];
bar.barConnect_Mat=[bar.barConnect_Mat; 12,13];
bar.barConnect_Mat=[bar.barConnect_Mat; 13,14];

bar.barConnect_Mat=[bar.barConnect_Mat; 15,16];
bar.barConnect_Mat=[bar.barConnect_Mat; 16,17];
bar.barConnect_Mat=[bar.barConnect_Mat; 17,18];

bar.barConnect_Mat=[bar.barConnect_Mat; 1,6];
bar.barConnect_Mat=[bar.barConnect_Mat; 10,15];

bar.E_Vec=[bar.E_Vec;Eeq*ones(16,1)];
bar.A_Vec=[bar.A_Vec;AbeamLong*ones(16,1)];
bar.L0_Vec=[bar.L0_Vec;L5/3*ones(16,1)];

%% Define horizontal bar elements
bar.barConnect_Mat=[bar.barConnect_Mat; 2,6];
bar.barConnect_Mat=[bar.barConnect_Mat; 3,7];
bar.barConnect_Mat=[bar.barConnect_Mat; 4,8];
bar.barConnect_Mat=[bar.barConnect_Mat; 5,9];

bar.barConnect_Mat=[bar.barConnect_Mat; 15,11];
bar.barConnect_Mat=[bar.barConnect_Mat; 16,12];
bar.barConnect_Mat=[bar.barConnect_Mat; 17,13];
bar.barConnect_Mat=[bar.barConnect_Mat; 18,14];

bar.E_Vec=[bar.E_Vec;Eeq*ones(8,1)];
bar.A_Vec=[bar.A_Vec;AbeamHor*ones(8,1)];
bar.L0_Vec=[bar.L0_Vec;L2*ones(8,1)];

%% Define horizontal rotational springs
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;1,2,6,3];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;6,3,7,8];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;3,4,8,5];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;8,5,9,23];

rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;10,15,11,12];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;15,16,12,17];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;12,17,13,14];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;17,18,14,24];

rotSpr.rotSprK_Vec=[rotSpr.rotSprK_Vec;SprBeam*ones(8,1)];
rotSpr.theta_StressFree_Vec=[rotSpr.theta_StressFree_Vec;pi*ones(8,1)];

%% Define diagonal bars
bar.barConnect_Mat=[bar.barConnect_Mat; 6,3];
bar.barConnect_Mat=[bar.barConnect_Mat; 3,8];
bar.barConnect_Mat=[bar.barConnect_Mat; 8,5];

bar.barConnect_Mat=[bar.barConnect_Mat; 15,12];
bar.barConnect_Mat=[bar.barConnect_Mat; 12,17];
bar.barConnect_Mat=[bar.barConnect_Mat; 17,14];

bar.E_Vec=[bar.E_Vec;Eeq*ones(6,1)];
bar.A_Vec=[bar.A_Vec;AbeamDiag*ones(6,1)];
bar.L0_Vec=[bar.L0_Vec;sqrt(L2*L2+L5*L5/9)*ones(6,1)];

%% Define diagonal rotational springs
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;2,6,3,7];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;7,3,8,4];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;4,8,5,9];

rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;11,15,12,16];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;16,12,17,13];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;13,17,14,18];

rotSpr.rotSprK_Vec=[rotSpr.rotSprK_Vec;SprBeamDiag*ones(6,1)];
rotSpr.theta_StressFree_Vec=[rotSpr.theta_StressFree_Vec;pi*ones(6,1)];

%% Define bar for the panels (perimeter)
bar.barConnect_Mat=[bar.barConnect_Mat; 9,18];
bar.barConnect_Mat=[bar.barConnect_Mat; 23,24];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;L3*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 5,19];
bar.barConnect_Mat=[bar.barConnect_Mat; 14,20];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;L1*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 19,22];
bar.barConnect_Mat=[bar.barConnect_Mat; 20,25];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;L4*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 22,23];
bar.barConnect_Mat=[bar.barConnect_Mat; 24,25];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;(L1+L2)*ones(2,1)];

%% Define bar for inside the panel
bar.barConnect_Mat=[bar.barConnect_Mat; 5,22];
bar.barConnect_Mat=[bar.barConnect_Mat; 14,25];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;sqrt(L1*L1+L4*L4)*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 5,23];
bar.barConnect_Mat=[bar.barConnect_Mat; 14,24];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;sqrt(L2*L2+L4*L4)*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 9,23];
bar.barConnect_Mat=[bar.barConnect_Mat; 18,24];

bar.E_Vec=[bar.E_Vec;Eeq*ones(2,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(2,1)];
bar.L0_Vec=[bar.L0_Vec;L4*ones(2,1)];

bar.barConnect_Mat=[bar.barConnect_Mat; 9,21];
bar.barConnect_Mat=[bar.barConnect_Mat; 18,21];
bar.barConnect_Mat=[bar.barConnect_Mat; 21,24];
bar.barConnect_Mat=[bar.barConnect_Mat; 21,23];

bar.E_Vec=[bar.E_Vec;Eeq*ones(4,1)];
bar.A_Vec=[bar.A_Vec;Apanel*ones(4,1)];
bar.L0_Vec=[bar.L0_Vec;sqrt(L4*L4/4+L3*L3/4)*ones(4,1)];

%% Define diagonal rotational springs for panel
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;19,5,22,23];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;24,14,25,20];

rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;22,5,23,9];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;18,24,14,25];

rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;5,9,23,21];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;21,18,24,14];

rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;23,9,21,18];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;9,18,21,24];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;23,21,24,18];
rotSpr.rotSprIJKL_Mat=[rotSpr.rotSprIJKL_Mat;9,21,23,24];

rotSpr.rotSprK_Vec=[rotSpr.rotSprK_Vec;SprPanel*ones(10,1)];
rotSpr.theta_StressFree_Vec=[rotSpr.theta_StressFree_Vec;pi*ones(10,1)];



%% Initialize and Plot for investigation

assembly=Assembly_MEMS();
assembly.bar=bar;
assembly.node=node;
assembly.rotSpr=rotSpr;

assembly.InitializeAssembly();

plots=Plot_MEMS();
plots.displayRange=2*10^-3;
plots.displayRangeRatio=0.5;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber()
plots.Plot_Shape_BarNumber()
plots.Plot_Shape_SprNumber()


%% Assign Mass Properties
% density of materials for the beam structure
rhoPZT=8000;
rhoSiO2=2600;
rhoPt=21000;
rhoAl=2700;

% density of materials for panels
rhoSi=2300;
tSi=25*10^-6;

% Solve the mass
node.mass_Vec=zeros(25,1);

% Mass node for beams
totalBeamMass=L5*L2*(tSiO2*rhoSiO2+tPtTop*rhoPt...
    +tPtTop*rhoPt+tPZT*rhoPZT);
node.mass_Vec(1:18)=totalBeamMass/8;

% Mass node for panels
totalPanelMass=(L1*2+L2*2+L3)*L4*(rhoSi*tSi);
node.mass_Vec([5,9,18,14,20,25,24,23,22,19])=...
    node.mass_Vec([5,9,18,14,20,25,24,23,22,19])...
    +totalPanelMass/30; % node at edge
node.mass_Vec(21)=totalPanelMass*2/3; % node at center


%% Define panels for plotting

plots.panelConnection{1}=[1,2,6];
plots.panelConnection{2}=[2,3,6];
plots.panelConnection{3}=[3,6,7];
plots.panelConnection{4}=[3,7,8];
plots.panelConnection{5}=[3,4,8];
plots.panelConnection{6}=[4,5,8];
plots.panelConnection{7}=[8,9,5];

plots.panelConnection{8}=[10,11,15];
plots.panelConnection{9}=[11,12,15];
plots.panelConnection{10}=[15,16,12];
plots.panelConnection{11}=[16,17,12];
plots.panelConnection{12}=[12,13,17];
plots.panelConnection{13}=[13,14,17];
plots.panelConnection{14}=[17,18,14];

plots.panelConnection{15}=[19,20,25,22];


%% Frequency analysis for the mode shape
% Obtain the mass matrix and stiffness matrix

[Fvec,Kmat]=assembly.SolveFK(assembly.node.currentU_Mat);
Mmat=node.FindMassMat;

supp=[1,1,1,1;
      2,1,1,1;
      6,1,1,1;
      10,1,1,1;
      11,1,1,1;
      15,1,1,1];

[Kadj,Fadj]=ModKforSupp(Kmat,supp,Fvec);
[Madj]=ModMforSupp(Mmat,supp);

[Umode,frequencySquared]=eig(Kadj,Madj);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));

freq1=sqrt(freq(1))/pi/2;
freq2=sqrt(freq(2))/pi/2;
freq3=sqrt(freq(3))/pi/2;

Umode1=Umode(:,index(1));
Usize=size(Umode1,1);

Umode1=reshape(Umode1,3,Usize/3)'/10000;

Umode2=Umode(:,index(2));
Umode2=reshape(Umode2,3,Usize/3)'/10000;

Umode3=Umode(:,index(3));
Umode3=reshape(Umode3,3,Usize/3)'/10000;

plots.Plot_DeformedShape(zeros(size(Umode1)),Umode1)
plots.Plot_DeformedShape(zeros(size(Umode1)),Umode2)
plots.Plot_DeformedShape(zeros(size(Umode1)),Umode3)



%% Statically Fold to an Angle

foldAngleRate=1.08;

sf=Solver_NR_Folding();
sf.assembly=assembly;
sf.supp=[1,1,1,1;
              2,1,1,1;
              6,1,1,1;
              10,1,1,1;
              11,1,1,1;
              15,1,1,1;];

activeCreaseLeft=[1,2,3,4];
activeCreaseRight=[5,6,7,8];

sf.increStep=20;
sf.tol=1*10^-8;

sf.targetRot=rotSpr.theta_Current_Vec;
sf.targetRot(activeCreaseLeft)=foldAngleRate*pi;
sf.targetRot(activeCreaseRight)=foldAngleRate*pi;

Uhis=sf.Solve();
plots.Plot_DeformedShape(zeros(size(Umode1)),squeeze(Uhis(end,:,:)))

%% Solve Frequency after folding 

[Fvec,Kmat]=assembly.SolveFK(assembly.node.currentU_Mat);
Mmat=node.FindMassMat;

supp=[1,1,1,1;
      2,1,1,1;
      6,1,1,1;
      10,1,1,1;
      11,1,1,1;
      15,1,1,1];

[Kadj,Fadj]=ModKforSupp(Kmat,supp,Fvec);
[Madj]=ModMforSupp(Mmat,supp);

[Umode,frequencySquared]=eig(Kadj,Madj);
frequencySquared=diag(frequencySquared);
[freq,index]=sort(abs(frequencySquared));


freq1Fold=sqrt(freq(1))/pi/2;
freq2Fold=sqrt(freq(2))/pi/2;
freq3Fold=sqrt(freq(3))/pi/2;

Umode1Fold=Umode(:,index(1));
Usize=size(Umode1Fold,1);
Umode1Fold=reshape(Umode1Fold,3,Usize/3)'/10000;

Umode2Fold=Umode(:,index(2));
Umode2Fold=reshape(Umode2Fold,3,Usize/3)'/10000;

Umode3Fold=Umode(:,index(3));
Umode3Fold=reshape(Umode3Fold,3,Usize/3)'/10000;

plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
    squeeze(Uhis(end,:,:))+Umode1Fold)
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
    squeeze(Uhis(end,:,:))+Umode2Fold)
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)),...
    squeeze(Uhis(end,:,:))+Umode3Fold)

%% Solve for transient dynamic responses
caa=Solver_CAA_Dynamics();
caa.assembly=assembly;
caa.supp=[1,1,1,1;
              2,1,1,1;
              6,1,1,1;
              10,1,1,1;
              11,1,1,1;
              15,1,1,1;];

caa.alpha=100;
caa.beta=0;

omega=freq1*2*pi;
dampingRatio=caa.alpha/2/omega+caa.beta*omega/2;
fprintf('damping ratio is %d \n', dampingRatio);

caa.dt=2*10^-7;
totalTime=0.8;
step=round(totalTime/caa.dt);

% External forces
caa.Fext=zeros(step,Usize/3,3); 

% Sine wave input
time=(1:step)*caa.dt;
caa.rotSprTargetAngle=pi*ones(step,24);


% Generate sweep function
y = chirp(time,1500,totalSweepTime,4000,'linear', 90);
% plot(time,y);

ExciatationAngle=4;

% Actuate one leg to generate twisting motion
for i=1:length(activeCreaseLeft)
    caa.rotSprTargetAngle(:,activeCreaseLeft(i))=...
        (foldAngleRate*pi+ExciatationAngle/180*pi*y)';
end

for i=1:length(activeCreaseRight)
    caa.rotSprTargetAngle(:,activeCreaseRight(i))=...
        (foldAngleRate*pi+0*ExciatationAngle/180*pi*y)';
end

Uhis=caa.Solve();

% plots.Plot_DeformedShape(zeros(size(Umode1)),squeeze(Uhis(end,:,:)))
% plots.Plot_DeformedHis(Uhis(1:100:end,:,:))

%% Plot the twisting angle
angleTwist=zeros(step,1);
for i=1:step
    node1=squeeze(Uhis(i,25,:))+node.coordinates_Mat(25,:)';
    node2=squeeze(Uhis(i,22,:))+node.coordinates_Mat(22,:)';
    vector=node1-node2;    
    angleTwist(i)=atan(vector(3)/vector(2))*180/pi;

end

figure
hold on
plot(time,angleTwist)

% fileName="Twist Angle.mat";
% save(fileName,"time","angleTwist")