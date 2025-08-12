clear all;
clc;
close all;

tic

%% Define the geometry of thick origami
L=50*10^(-3);
alpha=30/180*pi;
t=10*10^(-3);
gap=0*10^(-3);
faceThickness=1*10^(-3);

% Stiffness parameters of the structure
E=2*10^9; % Young's modulus
v=0.3; % Poisson's Ratio

sprStiff=0.1;
zlsprStiff=1000000;

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

zlspr.node_ij_mat=[zlspr.node_ij_mat; 97 116; 98 117];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 112 136; 114 138];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 121 145; 123 147];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 131 154; 132 155];




%% Connect panel to create big panels
zlspr.node_ij_mat=[zlspr.node_ij_mat; 8 31; 9 32; 11 34; 12 35];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 44 67; 45 68; 47 70; 48 71];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 37 61; 39 63; 40 64; 42 66];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 73 91; 75 93; 76 94; 78 96];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 103 127; 105 129; 106 130; 108 132];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 98 121; 99 122; 101 124; 102 125];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 44 67; 45 68; 47 70; 48 71];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 134 157; 135 158; 137 160; 138 161];

zlspr.k_vec=zlsprStiff*ones(88,1);



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
rot_spr_4N.mv_factor_vec=100000*ones(24,1);
rot_spr_4N.mv_vec=ones(24,1);
rot_spr_4N.mv_vec([1 3 5 7 9 11 12 15 17 19 21 23])=0;



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

%% Setup the loading controller
sf=Solver_NR_Folding_4N;
sf.assembly=assembly;
sf.supp=[1,1,1,1;
         2,1,1,1;
         3,1,1,1;];

sf.targetRot=rot_spr_4N.theta_current_vec;
originalRot=rot_spr_4N.theta_current_vec;

totalStep=150;
targetRate=0.9;

Uhis=zeros(totalStep,168,3);
rotHis=zeros(totalStep,24);


for k=1:totalStep
    
    for i=[1 3 5 7 9 11 12 15 17 19 21 23]
        sf.targetRot(i)=sf.targetRot(i)-targetRate/totalStep*pi; 
    end
    
    for i=[2 4 6 8 10 13 14 16 18 20 22 24]
        sf.targetRot(i)=sf.targetRot(i)+targetRate/totalStep*pi;
    end
    
    if k<5
        sf.increStep=10;
    else
        sf.increStep=2;
    end
    
    sf.tol=1*10^-7;
    sf.iterMax=50;
    
    Utemp=sf.Solve();
    Uhis(k,:,:)=Utemp(end,:,:);

    rotHis(k,:)=rot_spr_4N.theta_current_vec';

end

toc

%% Plot the folding animation and deformed shape
% Plot the folding history and kinematics of the system
plots.fileName='UniThick_60_Folding.gif';
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedHis(Uhis(1:10:end,:,:))



%% Fold angle history comparison
% Compare simulated folding angle with theoretical folding angle

% find the theoretical fold angle
theoryFold=zeros(totalStep,2);
for i=1:totalStep
    theoryFold(i,1)=targetRate*i/totalStep*pi;
    theoryFold(i,2)=2*(atan(tan(theoryFold(i,1)/2)/cos(pi/3)));
end

figure
hold on

% Plot the theoretical fold angle
plot(pi-theoryFold(:,1),pi+theoryFold(:,2),Color=[0.5,0.5,0.5])
plot(pi-theoryFold(:,1),pi-theoryFold(:,2),Color=[0.5,0.5,0.5])
plot(pi-theoryFold(:,1),pi+theoryFold(:,1),Color=[0.5,0.5,0.5])
plot(pi-theoryFold(:,1),pi-theoryFold(:,1),Color=[0.5,0.5,0.5])

% Plot simulation fold angle
for i=1:24
    plot(rotHis(:,1),rotHis(:,i),':',LineWidth=2)
end
