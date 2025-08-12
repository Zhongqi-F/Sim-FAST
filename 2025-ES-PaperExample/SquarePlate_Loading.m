clear all;
clc;
close all;

%% Define the geometry of thick origami
L=70*10^(-3);
alpha=30/180*pi;
t=30*10^(-3);
gap=0*10^(-3);


% Stiffness parameters of the structure
faceThickness=2*10^(-3);
E=4*10^9; % Young's modulus
v=0.3; % Poisson's Ratio

sprStiff=1;
sprMVfactor=1000;
zlsprStiffhinge=215806;
zlsprStifflatch=119999;


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

% modify the structure for a rectangular
node.coordinates_mat([19 22 43 46],1)=node.coordinates_mat([19 22 43 46],1)+L/2;
node.coordinates_mat([103 106 127 130],1)=node.coordinates_mat([103 106 127 130],1)+L/2;

node.coordinates_mat([39 42 63 66],1)=node.coordinates_mat([39 42 63 66],1)-L/2;
node.coordinates_mat([123 126 147 150],1)=node.coordinates_mat([123 126 147 150],1)-L/2;

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

hingeZLNum=size(zlspr.node_ij_mat,1);


%% Connect panel to create big panels
zlspr.node_ij_mat=[zlspr.node_ij_mat; 8 31; 9 32; 11 34; 12 35];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 44 67; 45 68; 47 70; 48 71];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 37 61; 39 63; 40 64; 42 66];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 73 91; 75 93; 76 94; 78 96];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 103 127; 105 129; 106 130; 108 132];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 98 121; 99 122; 101 124; 102 125];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 44 67; 45 68; 47 70; 48 71];
zlspr.node_ij_mat=[zlspr.node_ij_mat; 134 157; 135 158; 137 160; 138 161];

zlspr.k_vec=zlsprStiffhinge*ones(80,1);
zlspr.k_vec((hingeZLNum+1):end)=zlspr.k_vec((hingeZLNum+1):end)*1000;


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


%% Assume all creases are locked
fold_lock_number=[
    4 23 5 24;
    2 25 3 26;
    19 43 21 45;
    10 29 11 30;

    13 32 14 33;
    17 40 18 41;    
    34 58 36 60;
    49 68 50 69;
    
    53 76 54 77;
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
    143 166 144 167;
    145 164 146 165;
    
    ];

lockNum=size(fold_lock_number,1);

for i=1:lockNum
        zlspr.node_ij_mat=[zlspr.node_ij_mat;
            fold_lock_number(i,1:2);
            fold_lock_number(i,3:4)];
        zlspr.k_vec=[zlspr.k_vec;
            zlsprStifflatch;
            zlsprStifflatch];
end 



%% Plot for investigation
assembly.Initialize_Assembly();

plots=Plot_ThickOrigami();
plots.displayRange=L*5;
plots.displayRangeRatio=0.2;
plots.viewAngle1=-40;
plots.assembly=assembly;

plots.Plot_Shape_NodeNumber();
plots.Plot_Shape_SprNumber();
plots.Plot_Shape_CSTNumber();
plots.Plot_Shape_ZLsprNumber();



%% Setup the loading controller
nr=Solver_NR_Loading;
nr.assembly=assembly;

supportIndex=[38 39 63 62 122 123 147 146 20 19 43 44 104 103 127 128]';
nr.supp=[supportIndex(1:8),ones(8,1),ones(8,1),ones(8,1);
         supportIndex(9:end),zeros(8,1),ones(8,1),ones(8,1);];
nr.increStep=1;

loadIndex=[94 96 76 78 ]';
totalForce = 50;
nodalForce = totalForce/4;

step=100;
maximumLoad=0;
maximumStep=0;
failedCSTMember=0;
failedZLMember=0;

for i=1:step

    nr.load=[loadIndex,zeros(4,1),zeros(4,1),-nodalForce*ones(4,1)*i];
    Utemp=nr.Solve();
    
    Nzlspr=size(zlspr.node_ij_mat,1);
    zlsprF=zeros(Nzlspr,1);
    for j=1:Nzlspr
        index1=zlspr.node_ij_mat(j,1);
        index2=zlspr.node_ij_mat(j,2);
        node1=node.coordinates_mat(index1,:)+Utemp(end,index1,:);
        node2=node.coordinates_mat(index2,:)+Utemp(end,index2,:); 
        node1=squeeze(node1);
        node2=squeeze(node2);
        zlsprF(j)=abs(zlspr.k_vec(j)*norm(node1'-node2'));
    end

    [bar_strain_mat,l_Mat,x_reshape,trans_mat] = ...
        cst.Solve_Bar_Strain(squeeze(Utemp(end,:,:)),node.coordinates_mat);
    [cst_strain_mat] = cst.Solve_CST_Strain(bar_strain_mat,trans_mat);

    E_mat=zeros(3,3);
    E_mat(1,1)=E/(1-v*v);
    E_mat(2,2)=E/(1-v*v);
    E_mat(1,2)=v*E/(1-v*v);
    E_mat(2,1)=v*E/(1-v*v);
    E_mat(3,3)=E/(1+v)/2;

    cst_stress_mat=cst_strain_mat*E_mat;
    cstNum=size(cst_stress_mat,1);
    cstVonMises=zeros(cstNum,1);
    for j=1:cstNum
        sig1=cst_stress_mat(j,1);
        sig2=cst_stress_mat(j,2);
        tau=cst_stress_mat(j,3);

        cstVonMises(j)=sqrt(sig1^2-sig1*sig2+sig2^2+3*tau^2);
    end

    Uhis(i,:,:)=Utemp(end,:,:);

    if max(zlsprF(1:hingeZLNum))>1000 || max(zlsprF(81:end))>500 || max(cstVonMises)>40000000
        maximumLoad=totalForce*i;
        maximumStep=i;
        failedZLhingeMember = find(zlsprF(1:hingeZLNum)>1000);
        failedZLLatchMember = find(zlsprF(81:end)>500)+80;
        failedCSTMember = find(cstVonMises>40000000);
        maximumHingeForce=max(zlsprF(1:hingeZLNum));
        maximumLatchForce=max(zlsprF(81:end));
        maximumCSTStress=max(cstVonMises);
        break
    end    
end


plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
hold on
if isempty(failedZLLatchMember)
    nodeIndex1=zlspr.node_ij_mat(failedZLhingeMember,1);
    x=node.coordinates_mat(nodeIndex1,:)'+squeeze(Uhis(end,nodeIndex1,:));
else
    nodeIndex1=zlspr.node_ij_mat(failedZLLatchMember,1);
    x=node.coordinates_mat(nodeIndex1,:)'+squeeze(Uhis(end,nodeIndex1,:));
end

scatter3(x(1),x(2),x(3),500);
hold off


Urefhis=squeeze(Uhis(:,loadIndex,3));
Urefhis=abs(mean(Urefhis,2))*1000;
forceHis=totalForce*(1:maximumStep)';
Urefhis=[0;Urefhis];
forceHis=[0;forceHis];
figure
plot(Urefhis,forceHis)
xlabel('displacement (mm)')
ylabel('force (N)')

if isempty(failedCSTMember)
    disp('Connector failed')
else
    disp('3D printed triangle failed')
end

disp(maximumHingeForce);
disp(maximumLatchForce);
disp(maximumCSTStress);
disp(failedZLhingeMember);
disp(failedZLLatchMember);

