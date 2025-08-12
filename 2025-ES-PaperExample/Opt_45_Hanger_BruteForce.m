clear all;
clc;
close all;

%% Define the geometry and material property of thick origami
L=50*10^(-3);
t=15*10^(-3);
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



%% Define the nodal coordinates for thick origami
for j=1:3
    for i=1:5
        if mod(i,2)==1
            node.coordinates_mat(1+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+gap/2,0];
            node.coordinates_mat(2+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2/tan(pi/8),L*(j-1)+gap/2,0];
            node.coordinates_mat(3+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+L-gap/2/tan(pi/8),0];
            
            node.coordinates_mat(4+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+gap/2,t];
            node.coordinates_mat(5+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2/tan(pi/8),L*(j-1)+gap/2,t];
            node.coordinates_mat(6+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+L-gap/2/tan(pi/8),t];


            node.coordinates_mat(7+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2/tan(pi/8),L*(j-1)+L-gap/2,0];
            node.coordinates_mat(8+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+gap/2/tan(pi/8),0];
            node.coordinates_mat(9+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+L-gap/2,0];
            
            node.coordinates_mat(10+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2/tan(pi/8),L*(j-1)+L-gap/2,t];
            node.coordinates_mat(11+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+gap/2/tan(pi/8),t];
            node.coordinates_mat(12+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+L-gap/2,t];
        else

            node.coordinates_mat(1+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+L-gap/2,0];
            node.coordinates_mat(2+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2/tan(pi/8),L*(j-1)+L-gap/2,0];
            node.coordinates_mat(3+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+gap/2/tan(pi/8),0];
            
            node.coordinates_mat(4+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+L-gap/2,t];
            node.coordinates_mat(5+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2/tan(pi/8),L*(j-1)+L-gap/2,t];
            node.coordinates_mat(6+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2,L*(j-1)+gap/2/tan(pi/8),t];

            
            node.coordinates_mat(7+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+gap/2,0];
            node.coordinates_mat(8+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+L-gap/2/tan(pi/8),0];
            node.coordinates_mat(9+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2/tan(pi/8),L*(j-1)+gap/2,0];
            
            node.coordinates_mat(10+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+gap/2,t];
            node.coordinates_mat(11+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+L-gap/2,L*(j-1)+L-gap/2/tan(pi/8),t];
            node.coordinates_mat(12+12*(i-1)+60*(j-1),:)=...
                [L*(i-1)+gap/2/tan(pi/8),L*(j-1)+gap/2,t];

        end    
    end
end
node.coordinates_mat=node.coordinates_mat(7:174,:);


%% Initialize assembly
assembly=Assembly_ThickOrigami();
assembly.node=node;
assembly.rot_spr_4N=rot_spr_4N;
assembly.cst=cst;
assembly.zlspr=zlspr;



%% Set up the panels
for i=1:28
    assembly.Add_Triangle_Panel(6*(i-1)+1,6*(i-1)+2,6*(i-1)+3,...
        6*(i-1)+4,6*(i-1)+5,6*(i-1)+6,E,faceThickness,v);
end



%% Define the connectors
zlspr.node_ij_mat=[zlspr.node_ij_mat;5 12; 6 10];
zlspr.node_ij_mat=[zlspr.node_ij_mat;8 14; 9 15];
zlspr.node_ij_mat=[zlspr.node_ij_mat;16 22; 17 24];
zlspr.node_ij_mat=[zlspr.node_ij_mat;20 26; 21 25];

zlspr.node_ij_mat=[zlspr.node_ij_mat;29 36; 30 34];
zlspr.node_ij_mat=[zlspr.node_ij_mat;32 38; 33 39];
zlspr.node_ij_mat=[zlspr.node_ij_mat;40 46; 41 48];
zlspr.node_ij_mat=[zlspr.node_ij_mat;44 50; 45 49];

zlspr.node_ij_mat=[zlspr.node_ij_mat;59 65; 60 64];
zlspr.node_ij_mat=[zlspr.node_ij_mat;62 69; 63 67];
zlspr.node_ij_mat=[zlspr.node_ij_mat;71 77; 72 78];
zlspr.node_ij_mat=[zlspr.node_ij_mat;73 79; 74 81];

zlspr.node_ij_mat=[zlspr.node_ij_mat;83 89; 84 88];
zlspr.node_ij_mat=[zlspr.node_ij_mat;86 93; 87 91];
zlspr.node_ij_mat=[zlspr.node_ij_mat;95 101; 96 102];
zlspr.node_ij_mat=[zlspr.node_ij_mat;97 103; 98 105];

zlspr.node_ij_mat=[zlspr.node_ij_mat;107 113; 108 112];
zlspr.node_ij_mat=[zlspr.node_ij_mat;116 122; 117 121];

zlspr.node_ij_mat=[zlspr.node_ij_mat;125 132; 126 130];
zlspr.node_ij_mat=[zlspr.node_ij_mat;128 134; 129 135];
zlspr.node_ij_mat=[zlspr.node_ij_mat;136 142; 137 144];
zlspr.node_ij_mat=[zlspr.node_ij_mat;140 146; 141 145];

zlspr.node_ij_mat=[zlspr.node_ij_mat;149 156; 150 154];
zlspr.node_ij_mat=[zlspr.node_ij_mat;152 158; 153 159];
zlspr.node_ij_mat=[zlspr.node_ij_mat;160 166; 161 168];



%% Connect the big panel
zlspr.node_ij_mat=[zlspr.node_ij_mat;1 55; 3 56; 4 58; 6 59];
zlspr.node_ij_mat=[zlspr.node_ij_mat;25 79; 27 80; 28 82; 30 83];
zlspr.node_ij_mat=[zlspr.node_ij_mat;49 103; 51 104; 52 106; 54 107];
zlspr.node_ij_mat=[zlspr.node_ij_mat;61 115; 63 116; 64 118; 66 119];

zlspr.node_ij_mat=[zlspr.node_ij_mat;88 142; 90 143; 85 139; 87 140];
zlspr.node_ij_mat=[zlspr.node_ij_mat;112 166; 114 167; 109 163; 111 164];

zlspr.k_vec=zlsprStiff_hinge*ones(74,1);

%% Rotational Spring
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 4 12 10 11];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 61 69 67 68];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 124 132 130 131];


rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 18 22 24 23];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 75 79 81 80];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 138 142 144 143];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 28 34 36 35];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 85 91 93 92];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 148 154 156 155];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 42 46 48 47];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 99 105 103 104];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 162 168 166 167];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 58 64 65 66];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 70 77 78 76];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 82 88 89 90];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 94 101 102 100];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 106 112 113 114];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 115 121 122 123];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 127 134 135 133];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 139 145 146 147];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 151 158 159 157];

rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 7 14 15 13];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 19 25 26 27];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 31 38 39 37];
rot_spr_4N.node_ijkl_mat=[rot_spr_4N.node_ijkl_mat; 43 49 50 51];

rot_spr_4N.rot_spr_K_vec=sprStiff*ones(25,1);
rot_spr_4N.mv_factor_vec=sprMVfactor*ones(25,1);
rot_spr_4N.mv_vec=ones(25,1);
rot_spr_4N.mv_vec([2 5 8 11 18 19 20 21 22 23 24 25])=0;


%% Plot for investigation
assembly.Initialize_Assembly();

plots=Plot_ThickOrigami();
plots.displayRange=0.3;
plots.displayRangeRatio=0.5;
plots.viewAngle1=-40;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber();
plots.Plot_Shape_SprNumber();
plots.Plot_Shape_CSTNumber();
plots.Plot_Shape_ZLsprNumber();

%% Locking Pair for Optimization

fold_lock_number=[
    2 9 3 7;
    65 72 66 70;
    122 129 123 127;
    13 19 14 21;
    76 82 77 84;
    133 139 134 141;
    26 33 27 31;
    89 96 90 94;
    146 153 147 151;
    37 43 38 45;
    100 106 101 108;
    157 163 158 165;
    56 62 57 61;
    69 75 68 74;
    80 86 81 85;
    93 99 92 98;
    104 110 105 109;
    119 125 120 124;
    132 138 131 137;
    144 148 143 149;
    156 162 155 161;
    12 18 11 17;
    23 29 24 28;
    36 42 35 41;
    47 53 48 52];

 
%% Newton Raphson Loading for the Folding Motion
nr=Solver_NR_Loading;
nr.assembly=assembly;
nr.supp=[1,1,1,1;
         150,1,1,1;
         22,1,1,1;];
totalForce=0.1;
loadIndex=47;
nr.load=[47,0,0,-totalForce];
nr.increStep=5;
nr.tol=10^-8;
nr.iterMax=30;

Uhis=nr.Solve();
UstoreNonLock=squeeze(Uhis(end,:,:));
UhisStoreNonLock=Uhis;
Zmotion=Uhis(end,loadIndex,3);
averageDisp=abs(mean(Zmotion));
stiff_nonLock=totalForce*nr.increStep/averageDisp;


%% Optimization step
totalOptStep=25;
lockCreaseNum=3;
selectedCrease=zeros(lockCreaseNum,1);
previousStiff=stiff_nonLock;
currentStiff=0;
lockLine=zeros(lockCreaseNum,2);
UstoreLock=zeros(lockCreaseNum,168,3);

for j=1:lockCreaseNum*2
    zlspr.node_ij_mat(74+j,:)=[1,1];
    zlspr.k_vec(74+j)=zlsprStiff_latch;
end

% Find the optimized hinge placement using brute force
locks=1:25;
possibleCombo=nchoosek(locks,3);

comboNum=size(possibleCombo,1);
stiffRecord=zeros(comboNum,1);

for i=1:comboNum       

    zlspr.node_ij_mat(74+1,:)=fold_lock_number(possibleCombo(i,1),1:2);
    zlspr.node_ij_mat(74+2,:)=fold_lock_number(possibleCombo(i,1),3:4);  
    
    zlspr.node_ij_mat(74+3,:)=fold_lock_number(possibleCombo(i,2),1:2);
    zlspr.node_ij_mat(74+4,:)=fold_lock_number(possibleCombo(i,2),3:4);  
    
    zlspr.node_ij_mat(74+5,:)=fold_lock_number(possibleCombo(i,3),1:2);
    zlspr.node_ij_mat(74+6,:)=fold_lock_number(possibleCombo(i,3),3:4);  

    assembly.Initialize_Assembly
    
    Uhis=nr.Solve();   
    Zmotion=Uhis(end,loadIndex,3);
    averageDisp=abs(mean(Zmotion));    
    currentStiff=totalForce*nr.increStep/averageDisp;

    stiffRecord(i)=currentStiff;

end

[stif, index]=max(stiffRecord);

zlspr.node_ij_mat(74+1,:)=fold_lock_number(possibleCombo(index,1),1:2);
zlspr.node_ij_mat(74+2,:)=fold_lock_number(possibleCombo(index,1),3:4);  

zlspr.node_ij_mat(74+3,:)=fold_lock_number(possibleCombo(index,2),1:2);
zlspr.node_ij_mat(74+4,:)=fold_lock_number(possibleCombo(index,2),3:4);  

zlspr.node_ij_mat(74+5,:)=fold_lock_number(possibleCombo(index,3),1:2);
zlspr.node_ij_mat(74+6,:)=fold_lock_number(possibleCombo(index,3),3:4);  

nr.increStep=10;
UhisMax=nr.Solve();

%% Plot for comparison
plots.viewAngle1=30;
plots.viewAngle2=20;

lockLine=[fold_lock_number(possibleCombo(index,1),2:3);
        fold_lock_number(possibleCombo(index,2),2:3);
        fold_lock_number(possibleCombo(index,3),2:3);];

% Plot the no-lock deformed shape
plots.Plot_DeformedShape_OptLock(10*UstoreNonLock,...
    "Deformation without lock (10X deformation)",[ ]);
% plot the deformed shape with best lock positions
plots.Plot_DeformedShape_OptLock(10*squeeze(UhisMax(end,:,:)),...
    "Deformation with lock (10X deformation)",lockLine);

% plot the loading curve comparison
loadHis=(1:nr.increStep)*totalForce;
dispHisLock=squeeze(UhisMax(:,loadIndex,3));
dispHisLock=mean(dispHisLock,2);
dispHisLock=-[0;dispHisLock];
loadHis=[0,loadHis];

figure
hold on
plot(dispHisLock,loadHis)
legend('With 3 Locks (Brute force)')
xlabel('displacement (m)')
ylabel('force (N)')