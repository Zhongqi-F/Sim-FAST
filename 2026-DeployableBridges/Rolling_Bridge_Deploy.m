clear all; 
close all; 
clc;

%% Define the Rolling Bridge Geometry
% Define the nodes
H=1100/1000; % Height of the truss
HA=1100/1000; % Height of the active bar
L=1600/1000; % Length of unit span
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
kspr=3*barE*I/L;

node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=Std_Elements_Bars;
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
cst.v_vec=panelv*ones(cstNum,1);
cst.E_vec=panelE*ones(cstNum,1);
cst.t_vec=panelt*ones(cstNum,1);

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