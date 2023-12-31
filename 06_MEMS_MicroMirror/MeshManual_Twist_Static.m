%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Developer: Yi Zhu
% Advisor: Evgueni T. Filipov
%
% Acknowledgement: We would like to acknowledge the prior works from
% Ke Liu and Glaucio H. Paulino for establishing shared versions of
% nonrigid origami simulators. Their works paved the way for the new
% origami simulator presented in this package. 
%
% Reference:
% [1] Yi Zhu, Evgueni T. Filipov (2021). 'Sequentially Working Origami 
%     Multi-Physics Simulator (SWOMPS): A Versatile Implementation',
%     ASME IDETC-CIE Conference. DETC2021-68042. 
% [2] Y. Zhu, E. T. Filipov (2021). 'Rapid Multi-Physic Simulation for 
%     Electro-Thermal Origami Robotic Systems'  International Journal of 
%     Mechanical Sciences, 202-203, 106537.
% [3] Y. Zhu, E. T. Filipov (2020). 'A Bar and Hinge Model for Simulating 
%     Bistability in Origami Structures with Compliant Creases' Journal of 
%     Mechanisms and Robotics, 021110-1. 
% [4] Y. Zhu, E. T. Filipov (2019). 'An Efficient Numerical Approach for 
%     Simulating Contact in Origami Assemblages.' Proc. R. Soc. A, 475: 
%     20190366.       
% [5] Y. Zhu, E. T. Filipov (2019). 'Simulating compliant crease origami 
%     with a bar and hinge model.' IDETC/CIE 2019. 97119. 
% [6] K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear        
%     structural analysis of origami assemblages using the MERLIN2      
%     software.' Origami^7. 
% [7] K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid   
%     origami - An efficient computational approach.' Proc. R. Soc. A 473: 
%     20170348. 
% [8] K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to   
%     capture highly nonlinear behavior of non-rigid origami.'           
%     Proceedings of IASS Annual Symposium 2016. 
%
%%%%%%%%%%%%%%%%%%%%%%  Active Origami Simulator  %%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize the solver
clear all;clc;close all;
ori=OrigamiSolver;
ori.panelInnerBarStart=1;
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

%% Add node of the micro-mirror

ori.newNode(1,:)=[-L5/3,L3/2+L2,0];
ori.newNode(2,:)=[0,L3/2+L2,0];
ori.newNode(3,:)=[L5/3,L3/2+L2,0];
ori.newNode(4,:)=[L5*2/3,L3/2+L2,0];
ori.newNode(5,:)=[L5,L3/2+L2,0];

ori.newNode(6,:)=[0,L3/2,0];
ori.newNode(7,:)=[L5/3,L3/2,0];
ori.newNode(8,:)=[L5*2/3,L3/2,0];
ori.newNode(9,:)=[L5,L3/2,0];

ori.newNode(10,:)=[-L5/3,-L3/2-L2,0];
ori.newNode(11,:)=[0,-L3/2-L2,0];
ori.newNode(12,:)=[L5/3,-L3/2-L2,0];
ori.newNode(13,:)=[L5*2/3,-L3/2-L2,0];
ori.newNode(14,:)=[L5,-L3/2-L2,0];

ori.newNode(15,:)=[0,-L3/2,0];
ori.newNode(16,:)=[L5/3,-L3/2,0];
ori.newNode(17,:)=[L5*2/3,-L3/2,0];
ori.newNode(18,:)=[L5,-L3/2,0];

ori.newNode(19,:)=[L5,L3/2+L2+L1,0];
ori.newNode(20,:)=[L5,-L3/2-L2-L1,0];

ori.newNode(21,:)=[L5+L4/2,0,0];

ori.newNode(22,:)=[L5+L4,L3/2+L2+L1,0];
ori.newNode(23,:)=[L5+L4,L3/2,0];
ori.newNode(24,:)=[L5+L4,-L3/2,0];
ori.newNode(25,:)=[L5+L4,-L3/2-L2-L1,0];


%% solve for the area

tPtTop=0.1*10^-6;
tPZT=1*10^-6;
tPtBTM=0.15*10^-6;
tSiO2=1*10^-6;

tAl=1*10^-6;

ESiO2=66*10^9;
EPt=168*10^9;
EPZT=70*10^9;
EAl=70*10^9;

tBeamTotal=tSiO2+tPtTop+tPZT+tPtBTM;
EAbeam=L2*(tSiO2*ESiO2+tPtTop*EPt+tPZT*EPZT+tPtBTM*EPt);

% Assume material is homogenized with 100 GPa
AbeamLong=EAbeam/(100*10^9)/2;
AbeamDiag=3*AbeamLong;
AbeamHor=3*AbeamLong;
Apanel=50*AbeamLong;


% Bending stiffness
AreaMomentBeam=(tSiO2*L2*(tSiO2/2)*ESiO2)+...
    (tPtBTM*L2*(tPtBTM/2+tSiO2)*EPt)+...
    (tPZT*L2*(tPZT/2+tPtBTM+tSiO2)*EPZT)+...
    (tPtTop*L2*(tPtTop/2+tPZT+tPtBTM+tSiO2)*EPt);
c=AreaMomentBeam/EAbeam;

EIcomp=L2*(ESiO2*(1/12*tSiO2^3+tSiO2*(c-tSiO2/2)^2)+...
    EPt*(1/12*tPtBTM^3+tPtBTM*(c-tPtBTM/2-tSiO2)^2)+...
    EPZT*(1/12*tPZT^3+tPZT*(c-tPZT/2-tPtBTM-tSiO2)^2)+...
    EPt*(1/12*tPtTop^3+tPtTop*(c-tPtTop/2-tPZT-tPtBTM-tSiO2)^2));


tPary=1*10^-6;
EPary=2*10^9;
EAbeamPary=L2*(tSiO2*ESiO2+tPtTop*EPt+tPZT*EPZT+tPtBTM*EPt+tPary*EPary);
AreaMomentBeamPary=(tSiO2*L2*(tSiO2/2)*ESiO2)+...
    (tPtBTM*L2*(tPtBTM/2+tSiO2)*EPt)+...
    (tPZT*L2*(tPZT/2+tPtBTM+tSiO2)*EPZT)+...
    (tPtTop*L2*(tPtTop/2+tPZT+tPtBTM+tSiO2)*EPt)+...
    (tPary*L2*(tPary/2+tPtTop+tPZT+tPtBTM+tSiO2)*EPary);
c=AreaMomentBeamPary/EAbeamPary;
EIcompPary=L2*(ESiO2*(1/12*tSiO2^3+tSiO2*(c-tSiO2/2)^2)+...
    EPt*(1/12*tPtBTM^3+tPtBTM*(c-tPtBTM/2-tSiO2)^2)+...
    EPZT*(1/12*tPZT^3+tPZT*(c-tPZT/2-tPtBTM-tSiO2)^2)+...
    EPt*(1/12*tPtTop^3+tPtTop*(c-tPtTop/2-tPZT-tPtBTM-tSiO2)^2)+...
    EPary*(1/12*tPary^3+tPary*(c-tPary/2-tPtTop-tPZT-tPtBTM-tSiO2)^2));




SprBeam=4*EIcomp/L5;
SprBeamDiag=3*SprBeam;
SprPanel=50*SprBeam;



%% Add the bar element
ori.AddBar(1,2,AbeamLong,L5/3);
ori.AddBar(2,3,AbeamLong,L5/3);
ori.AddBar(3,4,AbeamLong,L5/3);
ori.AddBar(4,5,AbeamLong,L5/3);

ori.AddBar(6,7,AbeamLong,L5/3);
ori.AddBar(7,8,AbeamLong,L5/3);
ori.AddBar(8,9,AbeamLong,L5/3);

ori.AddBar(10,11,AbeamLong,L5/3);
ori.AddBar(11,12,AbeamLong,L5/3);
ori.AddBar(12,13,AbeamLong,L5/3);
ori.AddBar(13,14,AbeamLong,L5/3);

ori.AddBar(15,16,AbeamLong,L5/3);
ori.AddBar(16,17,AbeamLong,L5/3);
ori.AddBar(17,18,AbeamLong,L5/3);

ori.AddBar(1,6,AbeamLong,L5/3);
ori.AddBar(10,15,AbeamLong,L5/3);

% Horizontal spring of beam
ori.AddHinge(1,2,6,3,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(6,3,7,8,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(3,4,8,5,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(8,5,9,23,AbeamHor,L2,SprBeam,pi);

ori.AddHinge(10,15,11,12,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(15,16,12,17,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(12,17,13,14,AbeamHor,L2,SprBeam,pi);
ori.AddHinge(17,18,14,24,AbeamHor,L2,SprBeam,pi);

% Diagonal spring of beam
ori.AddHinge(2,6,3,7,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);
ori.AddHinge(7,3,8,4,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);
ori.AddHinge(4,8,5,9,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);

ori.AddHinge(11,15,12,16,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);
ori.AddHinge(16,12,17,13,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);
ori.AddHinge(13,17,14,18,AbeamDiag,sqrt(L2*L2+L5*L5/9),SprBeamDiag,pi);

% Bar for panel
ori.AddBar(9,18,Apanel,L3);
ori.AddBar(23,24,Apanel,L3);

ori.AddBar(5,19,Apanel,L1);
ori.AddBar(14,20,Apanel,L1);

ori.AddBar(19,22,Apanel,L4);
ori.AddBar(20,25,Apanel,L4);

ori.AddBar(22,23,Apanel,L1+L2);
ori.AddBar(24,25,Apanel,L1+L2);


% Diagonal spring of Panel
ori.AddHinge(19,5,22,23,Apanel,sqrt(L1*L1+L4*L4),SprPanel,pi);
ori.AddHinge(22,5,23,9,Apanel,sqrt(L2*L2+L4*L4),SprPanel,pi);
ori.AddHinge(18,24,14,25,Apanel,sqrt(L2*L2+L4*L4),SprPanel,pi);
ori.AddHinge(24,14,25,20,Apanel,sqrt(L1*L1+L4*L4),SprPanel,pi);

ori.AddHinge(5,9,23,21,Apanel,L4,SprPanel,pi);
ori.AddHinge(21,18,24,14,Apanel,L4,SprPanel,pi);

ori.AddHinge(23,9,21,18,Apanel,sqrt(L4*L4/4+L3*L3/4),SprPanel,pi);
ori.AddHinge(9,18,21,24,Apanel,sqrt(L4*L4/4+L3*L3/4),SprPanel,pi);
ori.AddHinge(23,21,24,18,Apanel,sqrt(L4*L4/4+L3*L3/4),SprPanel,pi);
ori.AddHinge(9,21,23,24,Apanel,sqrt(L4*L4/4+L3*L3/4),SprPanel,pi);

% Plot the results for inspection
ori.viewAngle1=190;
ori.viewAngle2=15;
ori.displayRange=[-0.4*10^(-3),1.2*10^(-3),...
    -0.9*10^(-3),0.9*10^(-3),-0.2*10^(-3),1.8*10^(-3)]'; % plotting range

ori.plotBars=1;
ori.Plot_MeshedOrigami(); % Plot the unmeshed origami for inspection;

%% Other initialization
newNodeNum = size(ori.newNode);
newNodeNum = newNodeNum(1);
ori.currentAppliedForce = zeros(newNodeNum,3);   
ori.currentU = zeros(newNodeNum,3);
ori.currentSprZeroStrain = pi*logical(ori.sprK);
ori.continuingLoading=1;

ori.oldCreaseNum=24;
ori.oldCreaseType=2*ones(24,1);
ori.creaseRef=zeros(24,8);


%% Assign Mechanical Properties
ori.panelE=100*10^9; % Assume E=100 GPa
ori.creaseE=100*10^9; % Assume E=100 GPa
ori.panelPoisson=0.3;
ori.creasePoisson=0.3; 

% density of polyer for the beam structure
rhoPZT=8000;
rhoSiO2=2600;
rhoPt=21000;
rhoAl=2700;

% for panels
rhoSi=2300;
tSi=25*10^-6;
rho=1/5*(rhoPZT+rhoSiO2+rhoPt+rhoAl+rhoSi);

% Solve the mass
ori.nodalMass=zeros(25,1);

% Mass node for beams
totalBeamMass=L5*L2*(tSiO2*rhoSiO2+tPtTop*rhoPt+tPtTop*rhoPt+tPZT*rhoPZT);
ori.nodalMass(1:18)=totalBeamMass/8;


% Mass node for panels
totalPanelMass=(L1*2+L2*2+L3)*L4*(rhoSi*tSi);
ori.nodalMass([5,9,18,14,20,25,24,23,22,19])=ori.nodalMass([5,9,18,14,20,25,24,23,22,19])+totalPanelMass/30; % node at edge
ori.nodalMass(21)=totalPanelMass*2/3; % node at center



%% check stiffness matrix
StiffMat=ori.Solver_CalcK();


%% Frequency analysis
ori.densityCrease=rho;
ori.densityPanel=rho;

frequency=ControllerFrequencyAnalysis;
frequency.supp=[1,1,1,1;
              2,1,1,1;
              6,1,1,1;
              10,1,1,1;
              11,1,1,1;
              15,1,1,1;];
 
ori.Solver_Solve()

[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
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

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode1);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode2);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode3);



%% Statically Fold to an Angle
selfFold=ControllerSelfFolding();

foldAngleRate=1.08;

selfFold.supp=[1,1,1,1;
              2,1,1,1;
              6,1,1,1;
              10,1,1,1;
              11,1,1,1;
              15,1,1,1;];

activeCreaseLeft=[17,18,19,20];
activeCreaseRight=[21,22,23,24];

selfFold.increStep=20;
selfFold.tol=1*10^-5;

selfFold.targetRotZeroStrain=pi*ones(ori.oldCreaseNum,1);
selfFold.targetRotZeroStrain(activeCreaseLeft)=foldAngleRate*pi;
selfFold.targetRotZeroStrain(activeCreaseRight)=foldAngleRate*pi;

selfFold.videoOpen=0;
selfFold.plotOpen=1;

ori.compliantCreaseOpen=0;
ori.loadingController{1}={"SelfFold",selfFold};
ori.Solver_Solve();


%% Solve Frequency

[frequencySquared,Umode]=ori.Dynamic_FrequencyAnalysis(frequency);
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

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode1Fold);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode2Fold);

ori.Plot_DeformedShape(ori.newNode+ori.currentU,...
     ori.newNode+ori.currentU+Umode3Fold);

%% Solve for static loading behavior

NR=ControllerNRLoading();
NR.supp=[1,1,1,1;
          2,1,1,1;
          6,1,1,1;
          10,1,1,1;
          11,1,1,1;
          15,1,1,1;];

NR.load=[22,10^-6,0,-10^-6;
         25,-10^-6,0,10^-6;];


NR.increStep=150;
NR.tol=10^-7;
NR.iterMax=50;

NR.videoOpen=0;

ori.loadingController{1}={"NR",NR};
ori.Solver_Solve()
toc
 


%% Find folding angle
% panel fold
Uhis=NR.Uhis;

step=NR.increStep;
angle=zeros(step,1);
for i=1:step
    node1=squeeze(Uhis(i,25,:))+ori.newNode(25,:)';
    node2=squeeze(Uhis(i,22,:))+ori.newNode(22,:)';
    vector=node1-node2;    
    angle(i)=atan(vector(3)/vector(2))*180/pi;
    if vector(1)<0
        angle(i)=atan(vector(3)/vector(2))*180/pi+180;
    end
end

figure
hold on
plot(squeeze(Uhis(:,25,3)),(1:step)*10^-6);

figure
hold on
plot(angle,(1:step)*10^-6);