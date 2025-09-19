clear all; 
close all; 
clc;

%% Define the Rolling Bridge Geometry
% Define the nodes
H=1100/1000; % Height of the truss
HA=1100/1000; % Height of the active bar
L1=1600/1000; % Length of 1st and last units' span
L2=1600/1000; % Length of other units' span
l=530/1000;
W=1600/1000;

barA=0.0001; % cross section area (m^2)
barE=80*10^9; % Youngs' modulus (Pa)
panelE=10*10^6;  % Young's modulus 10MPa
panelv=0.3; % Poisson's Ratio 0.3
panelt=100/1000; % Thicknes 10cm

% Property of active truss
activeBarE=80*10^9;

N=8;

I=1/12*0.01^4;
kspr=3*barE*I/L2;

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
node.coordinates_mat=[];
for i=1:N
    if i==1
    node.coordinates_mat=[node.coordinates_mat;
        0       0,   0; % 1
        0       W,   0; % 2
        L1,     0,   0; % 3 
        L1,     W,   0; % 4    
        l,      0,   H; % 5
        L1-l,   0,   H; % 6
        l,      W,   H; % 7
        L1-l,   W,   H; % 8
        L1,     0,   HA; %active bar's node
        L1,     W,   HA; %active bar's node
        ];
    elseif i==N
    node.coordinates_mat=[node.coordinates_mat;
        L1+((N-2)*L2),        0,   HA; %active bar's node
        L1+((N-2)*L2),        W,   HA; %active bar's node
        L1+((N-2)*L2+L1),     0,   0; 
        L1+((N-2)*L2+L1),     W,   0;     
        l+((N-2)*L2+L1),      0,   H; 
        L1-l+((N-2)*L2+L1),   0,   H; 
        l+((N-2)*L2+L1),      W,   H; 
        L1-l+((N-2)*L2+L1),   W,   H; 
        ];
    end
    if i>=3 && i<=N-1
    node.coordinates_mat=[node.coordinates_mat;
        L2+L1+(i-3)*L2,      0,   HA; %active bar's node
        L2+L1+(i-3)*L2,      W,   HA; %active bar's node
    ];
    end
    if i>=2 && i<=N-1   
    node.coordinates_mat=[node.coordinates_mat;
        L2+L1+(i-2)*L2,      0,   0;
        L2+L1+(i-2)*L2,      W,   0;    
        l+L1+(i-2)*L2,       0,   H; 
        (L2-l)+L1+(i-2)*L2,  0,   H; 
        l+L1+(i-2)*L2,       W,   H; 
        (L2-l)+L1+(i-2)*L2,  W,   H; 
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
        1  2  3;
        2  3  4;
        ];
    elseif i==N
    cst.node_ijk_mat=[cst.node_ijk_mat
        51  52  59;
        52  59  60;
        ];
    else
    cst.node_ijk_mat=[cst.node_ijk_mat
        3+(i-2)*8  4+(i-2)*8   11+(i-2)*8;
        4+(i-2)*8  11+(i-2)*8  12+(i-2)*8;
        ];
    end
end

cstNum=size(cst.node_ijk_mat,1);
cst.v_vec=panelv*ones(cstNum,1);
cst.E_vec=panelE*ones(cstNum,1);
cst.t_vec=panelt*ones(cstNum,1);

plots.Plot_Shape_CST_Number();


%% Define how normal bars are connected
% First we design the normal bar
for i=1:N
    if i==1
        bar.node_ij_mat=[bar.node_ij_mat;        
            1, 2;
            1, 3;
            2, 4;
            3, 4;
            1, 5;
            2, 7;
            3, 6;
            4, 8;
            5, 6;
            7, 8;
            6, 9;
            8, 10;
            ];
    elseif i==N
        bar.node_ij_mat=[bar.node_ij_mat;
            51, 59;
            59, 60;
            60, 52;
            51, 61;
            59, 62;
            52, 63;
            60, 64;
            61, 62;
            63, 64;
            61, 57;
            63, 58;
            ];
    else
        bar.node_ij_mat=[bar.node_ij_mat;
            3+(i-2)*8,   11+(i-2)*8;
            11+(i-2)*8,  12+(i-2)*8;
            12+(i-2)*8,  4+(i-2)*8;
            3+(i-2)*8,   13+(i-2)*8;
            13+(i-2)*8,  14+(i-2)*8;
            14+(i-2)*8,  11+(i-2)*8;
            4+(i-2)*8,   15+(i-2)*8;
            15+(i-2)*8,  16+(i-2)*8;
            16+(i-2)*8,  12+(i-2)*8;
            9+(i-2)*8,   13+(i-2)*8;
            10+(i-2)*8,  15+(i-2)*8;
            14+(i-2)*8,  17+(i-2)*8;
            16+(i-2)*8,  18+(i-2)*8;
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

% % Initialize the entire assembly 
% assembly.Initialize_Assembly();


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
rho_steel = 7850;         % Density of steel (kg/m^3)
g = 9.81;                 % Gravitational acceleration (m/s^2)

% Bars (including normal bars and actuator bars)
% Bar cross-section area
A_bar = barA;             % m^2, should match bar.A_vec

% 1. Total length of all normal bar elements
L_bar_total = 0;
barNodeMat = bar.node_ij_mat;
coords = node.coordinates_mat;
for i = 1:size(barNodeMat,1)
    n1 = barNodeMat(i,1);
    n2 = barNodeMat(i,2);
    p1 = coords(n1,:);
    p2 = coords(n2,:);
    len = norm(p1 - p2);
    L_bar_total = L_bar_total + len;
end

% 2. Total length of all actuator bar elements
L_actbar_total = 0;
if isfield(actBar,'node_ij_mat')
    actBarNodeMat = actBar.node_ij_mat;
    for i = 1:size(actBarNodeMat,1)
        n1 = actBarNodeMat(i,1);
        n2 = actBarNodeMat(i,2);
        p1 = coords(n1,:);
        p2 = coords(n2,:);
        len = norm(p1 - p2);
        L_actbar_total = L_actbar_total + len;
    end
end

% 3. Total bar length (all bars)
L_total = L_bar_total + L_actbar_total;

% 4. Total weight of all bars (Newtons)
W_bar = A_bar * L_total * rho_steel * g;

% CST panels
A_cst_total = 0;
cstNodeMat = cst.node_ijk_mat;
for i = 1:size(cstNodeMat,1)
    n1 = cstNodeMat(i,1);
    n2 = cstNodeMat(i,2);
    n3 = cstNodeMat(i,3);
    p1 = coords(n1,:);
    p2 = coords(n2,:);
    p3 = coords(n3,:);
    % Heron's formula for triangle area
    a = norm(p1 - p2);
    b = norm(p2 - p3);
    c = norm(p3 - p1);
    s = (a + b + c) / 2;
    area = sqrt(s * (s - a) * (s - b) * (s - c));
    A_cst_total = A_cst_total + area;
end

t_cst = panelt;  % Thickness of CST panel (m), matches cst.t_vec
% Weight of CST panels (Newtons)
W_cst = A_cst_total * t_cst * rho_steel * g;

% ----- Total self-weight -----
W_total = W_bar + W_cst;

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.2f m\n', L_total);
fprintf('Total area of all CST panels: %.2f m^2\n', A_cst_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Total CST panel weight: %.2f N\n', W_cst);
fprintf('Total self-weight of the bridge: %.2f N\n', W_total);
fprintf('-----------------------------\n');


%% Set up the self actuation solver
ta=Solver_NR_TrussAction;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

[T,K]=assembly.Solve_FK(zeros(size(node.coordinates_mat)));
spy(K)

ta.assembly=assembly;
ta.supp=[nodeNumVec,zeros(nodeNum,1),ones(nodeNum,1),zeros(nodeNum,1)];
ta.supp(1,2:4)=ones(1,3);
ta.supp(2,2:4)=ones(1,3);
ta.supp(3,2:4)=ones(1,3);
ta.supp(4,2:4)=ones(1,3);

% Set up the total loading step
ta.increStep=400;
% Set up the maximum iteration
ta.iterMax=30;
% Set up the tolorence
ta.tol=10^-4;

dL=0.4;
ta.targetL0=actBar.L0_vec;
ta.targetL0=ta.targetL0+dL;

% Solve for the deformation history
Uhis=ta.Solve();

% Plot the deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)));

% Also plot the deformation history
plots.fileName="Rolling_Bridge_Deploy.gif";
plots.Plot_Deformed_His(Uhis(1:10:end,:,:));