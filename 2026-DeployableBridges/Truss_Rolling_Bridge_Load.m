clear all;
clc;
close all;

H=1;
HA=1;
W=1;
L=1;
l=0.3;

N=4;

barA=0.01;
barE=2*10^9;
panel_E=200*10^9;
panel_t=0.05;
panel_v=0.2;
activeBarE=80*10^9;

% I=1/12*0.01^4;
% kspr=3*barE*I/L2*100000;

node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=CD_Elements_Bars;
cst=Vec_Elements_CST;

assembly=Assembly_Truss_Rolling_Bridge;
assembly.node=node;
assembly.bar=bar;
assembly.actBar=actBar;
assembly.cst=cst;

%% Define the nodal coordinates
node.coordinates_mat=[];
node.coordinates_mat=[node.coordinates_mat;
    0    0   0;
    L/2  0   H;
    L    0   0;
    0    W   0;
    L    W   0;
    L/2  W   H;
    L    0   H;
    L    W   H;
    ];

for i=2:N-1
    node.coordinates_mat=[node.coordinates_mat;
        L*i      0   0;
        L*i-L/2  0   H;
        L*i      W   0;
        L*i-L/2  W   H;
        L*i      0   H;
        L*i      W   H;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
    L*N      0   0;
    L*N-L/2  0   H;
    L*N      W   0;
    L*N-L/2  W   H;
    ];

% Set up the plotting function for inspection
plots=Plot_Truss_Rolling_Bridge();
plots.assembly=assembly;

% We will plot for the Rolling Bridge
plots.displayRange=[-0.5;N+0.5;-0.5;1.5;-0.5;1.5]; 

plots.viewAngle1=20;
plots.viewAngle2=20;


% Plot the nodal coordinates for inspection
plots.Plot_Shape_Node_Number()


%% Define how panels are designed
cst.node_ijk_mat=[cst.node_ijk_mat;
    1  3  4;
    3  4  5;
    ];

for i=2:N
    cst.node_ijk_mat=[cst.node_ijk_mat;
        3+(i-2)*6  5+(i-2)*6  9+(i-2)*6;
        5+(i-2)*6  9+(i-2)*6  11+(i-2)*6;
        ];
end    

cstNum=size(cst.node_ijk_mat,1);
cst.v_vec=panel_v*ones(cstNum,1);
cst.E_vec=panel_E*ones(cstNum,1);
cst.t_vec=panel_t*ones(cstNum,1);

plots.Plot_Shape_CST_Number();


%% Define how normal bars are connected
% First we design the normal bar
bar.node_ij_mat=[bar.node_ij_mat;
    1, 2;
    1, 3;
    2, 3;
    4, 5;
    4, 6;
    5, 6;
    2, 7;
    6, 8;
    1, 4;
    3, 5;
    3, 4;
    ];

for i=2:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        3+(i-2)*6,  9+(i-2)*6;
        3+(i-2)*6,  10+(i-2)*6;
        9+(i-2)*6,  10+(i-2)*6;
        5+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  12+(i-2)*6;
        11+(i-2)*6, 12+(i-2)*6;
        7+(i-2)*6,  10+(i-2)*6;
        8+(i-2)*6,  12+(i-2)*6;
        10+(i-2)*6, 13+(i-2)*6;
        12+(i-2)*6, 14+(i-2)*6;
        9+(i-2)*6,  11+(i-2)*6;
        5+(i-2)*6,  9+(i-2)*6;
        ];
end

bar.node_ij_mat=[bar.node_ij_mat;
    3+(N-2)*6,  9+(N-2)*6;
    3+(N-2)*6,  10+(N-2)*6;
    9+(N-2)*6,  10+(N-2)*6;
    5+(N-2)*6,  11+(N-2)*6;
    5+(N-2)*6,  12+(N-2)*6;
    11+(N-2)*6, 12+(N-2)*6;
    7+(N-2)*6,  10+(N-2)*6;
    8+(N-2)*6,  12+(N-2)*6;
    5+(N-2)*6,  9+(N-2)*6;
    9+(N-2)*6,  11+(N-2)*6;
    ];

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();

%% Define how actuator bars are connected
% Next we design the active bars
for i=1:N-1
    actBar.node_ij_mat=[actBar.node_ij_mat;
        3+(i-1)*6, 7+(i-1)*6;
        5+(i-1)*6, 8+(i-1)*6;
        ];
end

actBarNum=size(actBar.node_ij_mat,1);
actBar.A_vec=barA*ones(actBarNum,1);
actBar.E_vec=activeBarE*ones(actBarNum,1);

plots.Plot_Shape_ActBar_Number();

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
nr.supp(4,2:4)=ones(1,3);
nr.supp(6*N-3,2:4)=ones(1,3); % 21  45
nr.supp(6*N-1,2:4)=ones(1,3); % 23  47

force=1000;
step=10;
nr.load=[3*N-3,0,0,-force/2; % 9  21
         3*N-1,0,0,-force/2; % 11  23
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

% Also plot the deformation history
plots.fileName="Rolling_Bridge_Load.gif";
plots.Plot_Deformed_His(Uhis);

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
Uaverage=-mean(squeeze(Uhis(end,[3*N-3,3*N-1],3)));
Kstiff=step*force/Uaverage;


% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec)*barFailureForce/maxBarForce;
plots.Plot_Shape_Bar_Stress(bar_stress)



% Find the relationship betweeen the bar internal forces and load
loadatfail=force*step*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Normal bars: %d\n', barNum);
fprintf('Actuator bars: %d\n', actBarNum);
fprintf('Total bars: %d\n', barNum + actBarNum);