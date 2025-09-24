clear all; 
close all; 
clc;

%% Define the Rolling Bridge Geometry
% Define the nodes
H=1; %1118/1000; % Height of the truss
HA=1; %1055/1000; % Height of the active bar
L=1; %1315/1000; % Length of 1st and last units' span
l=0.3; %634/1000;
W=1; %1600/1000;

barA=0.01;
barE=2*10^9;
panel_E=200*10^9;
panel_t=0.05;
panel_v=0.2;
activeBarE=80*10^9;

% barA=0.0001; % cross section area (m^2)
% barE=80*10^9; % Youngs' modulus (Pa)
% panelE=10*10^6;  % Young's modulus 10MPa
% panelv=0.3; % Poisson's Ratio 0.3
% panelt=100/1000; % Thicknes 10cm

N=4;

I=1/12*0.01^4;
kspr=3*barE*I/L*100000;

node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=CD_Elements_Bars;
rot_spr_3N=CD_Elements_RotSprings_3N;
cst=Vec_Elements_CST;

assembly=Assembly_Rolling_Bridge;
assembly.node=node;
assembly.bar=bar;
assembly.actBar=actBar;
assembly.rot_spr_3N=rot_spr_3N;
assembly.cst=cst;


%% Define the nodal coordinates
node.coordinates_mat=[
    0     0,   0; 
    0     W,   0;
];

for i=1:N  

    if i<N
        node.coordinates_mat=[node.coordinates_mat;
            (i-1)*L+L,     0,   0; 
            (i-1)*L+L,     W,   0;   
            (i-1)*L+l,      0,   H; 
            (i-1)*L+L-l,   0,   H; 
            (i-1)*L+l,      W,   H; 
            (i-1)*L+L-l,   W,   H; 
            (i-1)*L+L,     0,   HA; 
            (i-1)*L+L,     W,   HA; 
            ];
    else 
        node.coordinates_mat=[node.coordinates_mat;
            (N-1)*L+L,     0,   0; 
            (N-1)*L+L,     W,   0;   
            (N-1)*L+l,      0,   H; 
            (N-1)*L+L-l,   0,   H; 
            (N-1)*L+l,      W,   H; 
            (N-1)*L+L-l,   W,   H; 
            ];
    end

end

% Set up the plotting function for inspection
plots=Plot_Rolling_Bridge();
plots.assembly=assembly;

% We will plot for the Rolling Bridge
plots.displayRange=[-2;14;-1;2;-1;10]; 

% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Define how panels are designed
for i=1:N
    if i==1
        cst.node_ijk_mat=[cst.node_ijk_mat
            (i-1)*8+1  (i-1)*8+2  (i-1)*8+3;
            (i-1)*8+2  (i-1)*8+3  (i-1)*8+4;
            ];
    else 
        cst.node_ijk_mat=[cst.node_ijk_mat
            (i-2)*8+3  (i-2)*8+4  (i-2)*8+11;
            (i-2)*8+4  (i-2)*8+11  (i-2)*8+12;
            ];
    end
end


cstNum=size(cst.node_ijk_mat,1);
cst.v_vec=panel_v*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.t_vec=panel_t*ones(cstNum,1);

plots.Plot_Shape_CST_Number();


%% Define how normal bars are connected
% First we design the normal bar
for i=1:N
    if i==1
    bar.node_ij_mat=[
            1 3;
            3 6;
            6 5;
            5 1;
            2 4;
            4 8;
            8 7;
            7 2;

            6 9;
            8 10;

            1 2;
        ];
    elseif i==N
    bar.node_ij_mat=[bar.node_ij_mat;
        (i-2)*8+3,   (i-2)*8+11;
        (i-2)*8+11,   (i-2)*8+14;
        (i-2)*8+13,   (i-2)*8+14;
        (i-2)*8+13,   (i-2)*8+3;
        
        (i-2)*8+4,   (i-2)*8+12;
        (i-2)*8+12,   (i-2)*8+16;
        (i-2)*8+15,   (i-2)*8+16;
        (i-2)*8+15,   (i-2)*8+4;

        (i-2)*8+9,   (i-2)*8+13;
        (i-2)*8+10,   (i-2)*8+15;   

        (i-2)*8+3,   (i-2)*8+4;
        (i-2)*8+11,   (i-2)*8+12;
        ];
    else
    bar.node_ij_mat=[bar.node_ij_mat;
        (i-2)*8+3,   (i-2)*8+11;
        (i-2)*8+11,   (i-2)*8+14;
        (i-2)*8+13,   (i-2)*8+14;
        (i-2)*8+13,   (i-2)*8+3;
        
        (i-2)*8+4,   (i-2)*8+12;
        (i-2)*8+12,   (i-2)*8+16;
        (i-2)*8+15,   (i-2)*8+16;
        (i-2)*8+15,   (i-2)*8+4;

        (i-2)*8+9,   (i-2)*8+13;
        (i-2)*8+14,   (i-2)*8+17;
        (i-2)*8+10,   (i-2)*8+15;
        (i-2)*8+16,   (i-2)*8+18;        

        (i-2)*8+3,   (i-2)*8+4;
        ];
    end

end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();

%% Define how actuator bars are connected
% Next we design the active bars
for i=1:N-1
    actBar.node_ij_mat=[actBar.node_ij_mat;
       3+(i-1)*8, 9+(i-1)*8;
       4+(i-1)*8, 10+(i-1)*8; 
       ];
end

actBarNum=size(actBar.node_ij_mat,1);
actBar.A_vec=barA*ones(actBarNum,1);
actBar.E_vec=activeBarE*ones(actBarNum,1);

plots.Plot_Shape_ActBar_Number();

% Initialize the entire assembly 
assembly.Initialize_Assembly();


%% Define how rotational springs are connected
for i=1:N 
    if i==1
       rot_spr_3N.node_ijk_mat=[rot_spr_3N.node_ijk_mat;
           5   1   3;
           1   3   6; 
           3   6   5;
           6   5   1; 

           7   2   4;
           2   4   8; 
           4   8   7;
           8   7   2;
           ]; 
    else
       rot_spr_3N.node_ijk_mat=[rot_spr_3N.node_ijk_mat;
           13+(i-2)*8   3+(i-2)*8    11+(i-2)*8;
           3+(i-2)*8    11+(i-2)*8   14+(i-2)*8; 
           11+(i-2)*8   14+(i-2)*8   13+(i-2)*8;
           14+(i-2)*8   13+(i-2)*8   3+(i-2)*8; 

           15+(i-2)*8   4+(i-2)*8    12+(i-2)*8;
           4+(i-2)*8    12+(i-2)*8   16+(i-2)*8; 
           12+(i-2)*8   16+(i-2)*8   15+(i-2)*8;
           16+(i-2)*8   15+(i-2)*8   4+(i-2)*8];
    end
end

% Define rotational stiffness
rotSprNum=size(rot_spr_3N.node_ijk_mat,1);
rot_spr_3N.rot_spr_K_vec=kspr*ones(rotSprNum,1);

% Plot the rotational spring number
plots.Plot_Shape_Spr_Number();

% Initialize the entire assembly 
assembly.Initialize_Assembly();


%% Calculate Self-weight of the Bridge
% Steel properties
rho_steel=7850;         % Density of steel (kg/m^3)
g=9.81;                 % Gravitational acceleration (m/s^2)

% Bars (including normal bars and actuator bars)
% Bar cross-section area
A_bar=barA;             % m^2, should match bar.A_vec

% 1. Total length of all normal bar elements
L_bar_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;
for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);
    L_bar_total=L_bar_total+len;
end

% 2. Total length of all actuator bar elements
L_actbar_total=0;
if isfield(actBar,'node_ij_mat')
    actBarNodeMat=actBar.node_ij_mat;
    for i=1:size(actBarNodeMat,1)
        n1=actBarNodeMat(i,1);
        n2=actBarNodeMat(i,2);
        p1=coords(n1,:);
        p2=coords(n2,:);
        len=norm(p1-p2);
        L_actbar_total=L_actbar_total+len;
    end
end

% 3. Total bar length (all bars)
L_total=L_bar_total+L_actbar_total;

% 4. Total weight of all bars (Newtons)
W_bar=A_bar*L_total*rho_steel*g;

% CST panels
A_cst_total=0;
cstNodeMat=cst.node_ijk_mat;
for i=1:size(cstNodeMat,1)
    n1=cstNodeMat(i,1);
    n2=cstNodeMat(i,2);
    n3=cstNodeMat(i,3);
    p1=coords(n1,:);
    p2=coords(n2,:);
    p3=coords(n3,:);
    % Heron's formula for triangle area
    a=norm(p1-p2);
    b=norm(p2-p3);
    c=norm(p3-p1);
    s=(a+b+c)/2;
    area=sqrt(s*(s-a)*(s-b)*(s-c));
    A_cst_total=A_cst_total+area;
end

t_cst=panel_t;  % Thickness of CST panel (m), matches cst.t_vec
% Weight of CST panels (Newtons)
W_cst=A_cst_total*t_cst*rho_steel*g;

% ----- Total self-weight -----
W_total=W_bar+W_cst;

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.2f m\n', L_total);
fprintf('Total area of all CST panels: %.2f m^2\n', A_cst_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Total CST panel weight: %.2f N\n', W_cst);
fprintf('Total self-weight of the bridge: %.2f N\n', W_total);
fprintf('-----------------------------\n');


%% Set up the self actuation solver
nr=Solver_NR_Loading;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

[T,K]=assembly.Solve_FK(zeros(size(node.coordinates_mat)));
spy(K)

nr.assembly=assembly;
nr.supp=[nodeNumVec,zeros(nodeNum,1),ones(nodeNum,1),zeros(nodeNum,1)];
nr.supp(1,2:4)=ones(1,3);
nr.supp(2,2:4)=ones(1,3);
nr.supp(8*N-5,2:4)=ones(1,3);
nr.supp(8*N-4,2:4)=ones(1,3);

force=1000;
step=10;
nr.load=[8*N/2-5,0,0,-force/2;
         8*N/2-4,0,0,-force/2;
         ];

% Set up the total loading step
nr.increStep=step;
% Set up the maximum iteration
nr.iterMax=30;
% Set up the tolorence
nr.tol=10^-5;


% Solve for the deformation history
Uhis=nr.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(20*squeeze(Uhis(end,:,:)));


% Failure load computation
truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);


% Find the maximum bar internal force
[maxBarForce,index]=max(abs(internal_force));

% Find failure force for the bar
sigma_u=300*10^6;
barFailureForce=sigma_u*barA;

% Find total bar length
barLtotal=sum(bar.L0_vec);

% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[11,12],3)));
Kstiff=step*force/Uaverage;

% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec)*barFailureForce/maxBarForce;
plots.Plot_Shape_Bar_Stress(bar_stress)

% Find the relationship betweeen the bar internal forces and load
loadatfail=step*force*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Normal bars: %d\n', barNum);
fprintf('Actuator bars: %d\n', actBarNum);
fprintf('Total bars: %d\n', barNum + actBarNum);