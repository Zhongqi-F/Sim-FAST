clear all;
clc;
close all;

%% Define the geometry of thick origami
L=150*10^(-3);
alpha=30/180*pi;
t=12*10^(-3);
gap=0*10^(-3);


% Stiffness parameters of the structure
faceThickness=1*10^(-3);
% Please note that this sample has a smaller wall thickness compared to the
% plate used in the four-point bending test. This is because the plate is
% also thinner.

E=4*10^9; % Young's modulus
v=0.3; % Poisson's Ratio

sprStiff=1;
sprMVfactor=1000;
zlsprStiff_latch=119999;
zlsprStiff_hinge=215806;


% Initialize Elements
node=Elements_Nodes;
rot_spr_4N=Vec_Elements_RotSprings_4N_Directional;
cst=Vec_Elements_CST;
zlspr=Vec_Elements_Zero_L_Spring;



%% thick origami geometry
for i=1:3
    node.coordinates_mat(1+6*(i-1),:)=...
        [gap*sin(pi/3)+(L+2*gap*sin(pi/3))*(i-1),0,0];
    node.coordinates_mat(2+6*(i-1),:)=...
        [gap*sin(pi/3)+L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3),0];
    node.coordinates_mat(3+6*(i-1),:)=...
        [gap*sin(pi/3)+L+(L+2*gap*sin(pi/3))*(i-1),0,0];
    
    node.coordinates_mat(4+6*(i-1),:)=...
        [gap*sin(pi/3)+(L+2*gap*sin(pi/3))*(i-1),0,t];
    node.coordinates_mat(5+6*(i-1),:)=...
        [gap*sin(pi/3)+L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3),t];
    node.coordinates_mat(6+6*(i-1),:)=...
        [gap*sin(pi/3)+L+(L+2*gap*sin(pi/3))*(i-1),0,t];
end

for i=1:4
    node.coordinates_mat(18+1+6*(i-1),:)=...
        [-L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2,0];
    node.coordinates_mat(18+2+6*(i-1),:)=...
        [(L+2*gap*sin(pi/3))*(i-1),gap/2,0];
    node.coordinates_mat(18+3+6*(i-1),:)=...
        [L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2,0];
    
    node.coordinates_mat(18+4+6*(i-1),:)=...
        [-L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2,t];
    node.coordinates_mat(18+5+6*(i-1),:)=...
        [(L+2*gap*sin(pi/3))*(i-1),gap/2,t];
    node.coordinates_mat(18+6+6*(i-1),:)=...
        [L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2,t];
end

for i=1:4
    node.coordinates_mat(42+1+6*(i-1),:)=...
        [-L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2*3,0];
    node.coordinates_mat(42+2+6*(i-1),:)=...
        [(L+2*gap*sin(pi/3))*(i-1),gap/2*3+2*L*sin(pi/3),0];
    node.coordinates_mat(42+3+6*(i-1),:)=...
        [L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2*3,0];
    
    node.coordinates_mat(42+4+6*(i-1),:)=...
        [-L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2*3,t];
    node.coordinates_mat(42+5+6*(i-1),:)=...
        [(L+2*gap*sin(pi/3))*(i-1),gap/2*3+2*L*sin(pi/3),t];
    node.coordinates_mat(42+6+6*(i-1),:)=...
        [L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+gap/2*3,t];
end

for i=1:3
    node.coordinates_mat(66+1+6*(i-1),:)=...
        [gap*sin(pi/3)+(L+2*gap*sin(pi/3))*(i-1),2*L*sin(pi/3)+2*gap,0];
    node.coordinates_mat(66+2+6*(i-1),:)=...
        [gap*sin(pi/3)+L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+2*gap,0];
    node.coordinates_mat(66+3+6*(i-1),:)=...
        [gap*sin(pi/3)+L+(L+2*gap*sin(pi/3))*(i-1),2*L*sin(pi/3)+2*gap,0];
    
    node.coordinates_mat(66+4+6*(i-1),:)=...
        [gap*sin(pi/3)+(L+2*gap*sin(pi/3))*(i-1),2*L*sin(pi/3)+2*gap,t];
    node.coordinates_mat(66+5+6*(i-1),:)=...
        [gap*sin(pi/3)+L/2+(L+2*gap*sin(pi/3))*(i-1),L*sin(pi/3)+2*gap,t];
    node.coordinates_mat(66+6+6*(i-1),:)=...
        [gap*sin(pi/3)+L+(L+2*gap*sin(pi/3))*(i-1),2*L*sin(pi/3)+2*gap,t];
end

node.coordinates_mat(85:84+84,:)=node.coordinates_mat(1:84,:);
node.coordinates_mat(85:84+84,2)=node.coordinates_mat(85:84+84,2)...
    +2*L*sin(pi/3)+3*gap;


%% Initialize assembly
assembly=Assembly_ThickOrigami();
assembly.node=node;
assembly.rot_spr_4N=rot_spr_4N;
assembly.cst=cst;
assembly.zlspr=zlspr;


%% Set up the panels
for i=1:14*2
    assembly.Add_Triangle_Panel(6*(i-1)+1,6*(i-1)+2,6*(i-1)+3,...
        6*(i-1)+4,6*(i-1)+5,6*(i-1)+6,E,faceThickness,v);
end

%% Define the connectors
zlspr.node_ij_mat=[zlspr.node_ij_mat;1 20; 2 21];
zlspr.node_ij_mat=[zlspr.node_ij_mat;5 28; 6 29];
zlspr.node_ij_mat=[zlspr.node_ij_mat;22 46; 24 48];
zlspr.node_ij_mat=[zlspr.node_ij_mat;7 26; 8 27];

zlspr.node_ij_mat=[zlspr.node_ij_mat;16 35; 17 36];
zlspr.node_ij_mat=[zlspr.node_ij_mat;14 37; 15 38];
zlspr.node_ij_mat=[zlspr.node_ij_mat;31 55; 33 57];
zlspr.node_ij_mat=[zlspr.node_ij_mat;52 71; 53 72];

zlspr.node_ij_mat=[zlspr.node_ij_mat;50 73; 51 74];
zlspr.node_ij_mat=[zlspr.node_ij_mat;58 77; 59 78];
zlspr.node_ij_mat=[zlspr.node_ij_mat;61 80; 62 81];
zlspr.node_ij_mat=[zlspr.node_ij_mat;67 85; 69 87];

zlspr.node_ij_mat=[zlspr.node_ij_mat; 82 100; 84 102];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 88 107; 89 108];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 91 110; 92 111];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 95 118; 96 119];

zlspr.node_ij_mat=[zlspr.node_ij_mat; 97 116; 98 117];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 112 136; 114 138];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 121 145; 123 147];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 131 154; 132 155];

zlspr.node_ij_mat=[zlspr.node_ij_mat; 133 152; 134 153];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 142 161; 143 162];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 140 163; 141 164];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 148 167; 149 168];



%% Connect panel to create big panels
zlspr.node_ij_mat=[zlspr.node_ij_mat; 8 31; 9 32; 11 34; 12 35];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 44 67; 45 68; 47 70; 48 71];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 37 61; 39 63; 40 64; 42 66];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 73 91; 75 93; 76 94; 78 96];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 103 127; 105 129; 106 130; 108 132];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 98 121; 99 122; 101 124; 102 125];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 134 157; 135 158; 137 160; 138 161];

zlspr.k_vec=zlsprStiff_hinge*ones(76,1);
zlspr.k_vec(49:76,1)=zlspr.k_vec(49:76,1)*10; % This is the quad panel



%% Set up the rotational hinges
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 3 20 21 19];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 4 28 29 30];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 9 26 27 25];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 18 35 36 34];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 13 37 38 39];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 23 46 48 47];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 32 55 57 56];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 54 71 72 70];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 49 73 74 75];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 60 77 78 76];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 63 80 81 79];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 68 85 87 86];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 83 100 102 101];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 90 107 108 106];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 93 110 111 109];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 94 118 119 120];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 99 116 117 115];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 113 136 138 137];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 122 145 147 146];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 130 131 132 156];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 151 152 153 135];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 160 161 162 144];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 139 140 141 165];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 166 167 168 150];

rot_spr_4N.rot_spr_K_vec=sprStiff*ones(24,1);
rot_spr_4N.mv_factor_vec=sprMVfactor*ones(24,1);
rot_spr_4N.mv_vec=ones(24,1);
rot_spr_4N.mv_vec([1 3 5 7 9 11 12 15 17 19 21 23])=0;


%% Plot for investigation
assembly.Initialize_Assembly();

plots=Plot_ThickOrigami();
plots.displayRange=1;
plots.displayRangeRatio=0.5;
plots.viewAngle1=-40;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber();
plots.Plot_Shape_SprNumber();
plots.Plot_Shape_CSTNumber();
plots.Plot_Shape_ZLsprNumber();

%% Locking Pair for Optimization
fold_lock_number=[
    4 23 5 24;
    2 25 3 26;
    10 29 11 30;
    13 32 14 33;
    18 41 17 40;
    19 43 21 45;
    34 58 36 60;
    49 68 50 69;
    54 77 53 76;
    55 74 56 75;
    64 83 65 84;
    70 88 72 90;
    79 97 81 99;
    85 104 86 105;
    94 113 95 114;
    92 115 93 116;
    100 119 101 120;
    109 133 111 135;
    124 148 126 150;
    129 152 128 151;
    136 155 137 156;
    139 158 140 159;
    144 167 143 166;
    145 164 146 165];


    
%% Setup the loading controller
nr=Solver_NR_Loading;
nr.assembly=assembly;    

supportIndex=[19 29 21 22 23 24]';
nr.supp=[supportIndex,ones(6,1),ones(6,1),ones(6,1)];

loadIndex=[148; 149; 150];
totalForce = 0.05;
nodalForce = totalForce/3;
nr.load=[loadIndex,zeros(3,1),zeros(3,1),-nodalForce*ones(3,1)];
nr.increStep= 5;
nr.tol=10^-6;
nr.iterMax=30;

Uhis=nr.Solve();
UstoreNonLock=squeeze(Uhis(end,:,:));
UhisStoreNonLock=Uhis;
Zmotion=Uhis(end,loadIndex,3);
averageDisp=abs(mean(Zmotion));
stiff_nonLock=totalForce*nr.increStep/averageDisp;

loadHisNonLock=(1:nr.increStep)*totalForce;

%% Optimization step
totalOptStep=24;
lockCreaseNum=5;
selectedCrease=zeros(lockCreaseNum,1);
previousStiff=stiff_nonLock;
currentStiff=0;
lockLine=zeros(lockCreaseNum,2);
UstoreLock=zeros(lockCreaseNum,168,3);

for j=1:lockCreaseNum*2
    zlspr.node_ij_mat(76+j,:)=[1,1];
    zlspr.k_vec(76+j)=zlsprStiff_latch;
end

% Find the optimized hinge placement
for i=1:lockCreaseNum
    
    locks=1:24;
    locks=locks';
    locks(ismember(locks,selectedCrease))=[];
    unusedLocks=locks;

    nodalForce = i*totalForce/3;
    nr.load=[loadIndex,zeros(3,1),zeros(3,1),-nodalForce*ones(3,1)];

    for j=unusedLocks'

        previousLock=zlspr.node_ij_mat((76+2*i-1):(76+2*i),:);
            
        zlspr.node_ij_mat(76+2*i-1,:)=fold_lock_number(j,1:2);
        zlspr.node_ij_mat(76+2*i,:)=fold_lock_number(j,3:4);  

        assembly.Initialize_Assembly
        
        Uhis=nr.Solve();   
        Zmotion=Uhis(end,loadIndex,3);
        averageDisp=abs(mean(Zmotion));    
        currentStiff=i*totalForce*nr.increStep/averageDisp;

        if currentStiff<previousStiff
            zlspr.node_ij_mat((76+2*i-1):(76+2*i),:)=previousLock;
        else
            previousStiff=currentStiff;   
            selectedCrease(i)=j;
            lockLine(i,:)=fold_lock_number(j,2:3);
            UstoreLock(i,:,:)=squeeze(Uhis(end,:,:));
            UhisStoreLock{i}=Uhis;
        end
    end
end

%% Find a random lock case for comparison
randomLock=randperm(24,lockCreaseNum);

nodalForce = 2*totalForce/3;
nr.load=[loadIndex,zeros(3,1),zeros(3,1),-nodalForce*ones(3,1)];

for i=1:lockCreaseNum
    zlspr.node_ij_mat(76+2*i-1,:)=fold_lock_number(randomLock(i),1:2);
    zlspr.node_ij_mat(76+2*i,:)=fold_lock_number(randomLock(i),3:4);  
end

assembly.Initialize_Assembly
Uhis=nr.Solve();    
dispHisRandomLock=squeeze(Uhis(:,loadIndex,3));
dispHisRandomLock=mean(dispHisRandomLock,2);

loadHisRandom=(1:nr.increStep)*totalForce*2;

%% Plot for comparison

plots.viewAngle1=0;
plots.viewAngle2=20;

% Plot the no-lock deformed shape
plots.Plot_DeformedShape_OptLock(10*UstoreNonLock,...
    "Deformation without lock (10X deformation)",[ ]);
% plot the deformed shape with applied locks
plots.Plot_DeformedShape_OptLock(10*squeeze(UstoreLock(end,:,:)),...
    "Deformation with lock (10X deformation)",lockLine);

% plot the loading curve comparison
dispHisNonLock=squeeze(UhisStoreNonLock(:,loadIndex,3));
dispHisNonLock=mean(dispHisNonLock,2);

UhisStoreLock5=UhisStoreLock{5};
dispHisLock5=squeeze(UhisStoreLock5(:,loadIndex,3));
dispHisLock5=mean(dispHisLock5,2);

UhisStoreLock3=UhisStoreLock{3};
dispHisLock3=squeeze(UhisStoreLock3(:,loadIndex,3));
dispHisLock3=mean(dispHisLock3,2);

loadHis5=(1:nr.increStep)*totalForce*5;
loadHis3=(1:nr.increStep)*totalForce*3;

dispHisLock5=-[0;dispHisLock5]*1000;
dispHisLock3=-[0;dispHisLock3]*1000;
dispHisRandomLock=-[0;dispHisRandomLock]*1000;
dispHisNonLock=-[0;dispHisNonLock]*1000;

loadHis5=[0,loadHis5]';
loadHis3=[0,loadHis3]';
loadHisRandom=[0,loadHisRandom]';
loadHisNonLock=[0,loadHisNonLock]';

figure
hold on
plot(dispHisLock5,loadHis5)
plot(dispHisLock3,loadHis3)
plot(dispHisRandomLock,loadHisRandom)
plot(dispHisNonLock,loadHisNonLock)
legend('With 5 Locks (Optimized)','With 3 Locks (Optimized)'...
    ,'With 5 Locks (Random)','Without Locks')
xlabel('displacement (mm)')
ylabel('force (N)')
