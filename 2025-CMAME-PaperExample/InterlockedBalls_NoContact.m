clear all
close all
clc

tic

%% Define Geometry
R=0.1;
W=0.04;

node=Elements_Nodes;

%% First Sphere
% tip square
node.coordinates_mat=[node.coordinates_mat;
    -W/2    -W/2   R;
    -W/2    W/2    R;
    W/2     -W/2   R;
    W/2     W/2    R;];

node.coordinates_mat=[node.coordinates_mat;
    -W/2    -W/2   -R;
    -W/2    W/2    -R;
    W/2     -W/2   -R;
    W/2     W/2    -R;];

node.coordinates_mat=[node.coordinates_mat;
   -R   -W/2   -W/2;
   -R   -W/2   W/2;
   -R   W/2    -W/2;
   -R   W/2    W/2];

node.coordinates_mat=[node.coordinates_mat;
   R   -W/2   -W/2;
   R   -W/2   W/2;
   R   W/2    -W/2;
   R   W/2    W/2];

node.coordinates_mat=[node.coordinates_mat;
   -W/2   -R   -W/2;
   -W/2   -R   W/2;
   W/2    -R   -W/2;
   W/2    -R   W/2];

node.coordinates_mat=[node.coordinates_mat;
   -W/2   R   -W/2;
   -W/2   R   W/2;
   W/2    R   -W/2;
   W/2    R   W/2];

% first round
node.coordinates_mat=[node.coordinates_mat;
    -W/2    R*cos(1/3*pi/2)   R*sin(1/3*pi/2);
    W/2     R*cos(1/3*pi/2)   R*sin(1/3*pi/2);
    -W/2    R*cos(2/3*pi/2)   R*sin(2/3*pi/2);
    W/2     R*cos(2/3*pi/2)   R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    -W/2    R*cos(-1/3*pi/2)   R*sin(-1/3*pi/2);
    W/2     R*cos(-1/3*pi/2)   R*sin(-1/3*pi/2);
    -W/2    R*cos(-2/3*pi/2)   R*sin(-2/3*pi/2);
    W/2     R*cos(-2/3*pi/2)   R*sin(-2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    -W/2    -R*cos(1/3*pi/2)   -R*sin(1/3*pi/2);
    W/2     -R*cos(1/3*pi/2)   -R*sin(1/3*pi/2);
    -W/2    -R*cos(2/3*pi/2)   -R*sin(2/3*pi/2);
    W/2     -R*cos(2/3*pi/2)   -R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    -W/2    -R*cos(-1/3*pi/2)   -R*sin(-1/3*pi/2);
    W/2     -R*cos(-1/3*pi/2)   -R*sin(-1/3*pi/2);
    -W/2    -R*cos(-2/3*pi/2)   -R*sin(-2/3*pi/2);
    W/2     -R*cos(-2/3*pi/2)   -R*sin(-2/3*pi/2);];

% second round
node.coordinates_mat=[node.coordinates_mat;
    R*cos(1/3*pi/2)   -W/2    R*sin(1/3*pi/2);
    R*cos(1/3*pi/2)   W/2     R*sin(1/3*pi/2);
    R*cos(2/3*pi/2)   -W/2    R*sin(2/3*pi/2);
    R*cos(2/3*pi/2)   W/2     R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    R*cos(-1/3*pi/2)   -W/2    R*sin(-1/3*pi/2);
    R*cos(-1/3*pi/2)   W/2     R*sin(-1/3*pi/2);
    R*cos(-2/3*pi/2)   -W/2    R*sin(-2/3*pi/2);
    R*cos(-2/3*pi/2)   W/2     R*sin(-2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    -R*cos(1/3*pi/2)   -W/2    -R*sin(1/3*pi/2);
    -R*cos(1/3*pi/2)   W/2     -R*sin(1/3*pi/2);
    -R*cos(2/3*pi/2)   -W/2    -R*sin(2/3*pi/2);
    -R*cos(2/3*pi/2)   W/2     -R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    -R*cos(-1/3*pi/2)   -W/2    -R*sin(-1/3*pi/2);
    -R*cos(-1/3*pi/2)   W/2     -R*sin(-1/3*pi/2);
    -R*cos(-2/3*pi/2)   -W/2    -R*sin(-2/3*pi/2);
    -R*cos(-2/3*pi/2)   W/2     -R*sin(-2/3*pi/2);];



%% Second Sphere
% tip square
factor=6;
node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    -W/2   R;
    factor*R/4+-W/2    W/2    R;
    factor*R/4+W/2     -W/2   R;
    factor*R/4+W/2     W/2    R;];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    -W/2   -R;
    factor*R/4+-W/2    W/2    -R;
    factor*R/4+W/2     -W/2   -R;
    factor*R/4+W/2     W/2    -R;];

node.coordinates_mat=[node.coordinates_mat;
   factor*R/4+-R   -W/2   -W/2;
   factor*R/4+-R   -W/2   W/2;
   factor*R/4+-R   W/2    -W/2;
   factor*R/4+-R   W/2    W/2];

node.coordinates_mat=[node.coordinates_mat;
   factor*R/4+R   -W/2   -W/2;
   factor*R/4+R   -W/2   W/2;
   factor*R/4+R   W/2    -W/2;
   factor*R/4+R   W/2    W/2];

node.coordinates_mat=[node.coordinates_mat;
   factor*R/4+-W/2   -R   -W/2;
   factor*R/4+-W/2   -R   W/2;
   factor*R/4+W/2    -R   -W/2;
   factor*R/4+W/2    -R   W/2];

node.coordinates_mat=[node.coordinates_mat;
   factor*R/4+-W/2   R   -W/2;
   factor*R/4+-W/2   R   W/2;
   factor*R/4+W/2    R   -W/2;
   factor*R/4+W/2    R   W/2];

% first round
node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+R*cos(1/3*pi/2)    R*sin(1/3*pi/2)   -W/2;
    factor*R/4+R*cos(1/3*pi/2)     R*sin(1/3*pi/2)   W/2;
    factor*R/4+R*cos(2/3*pi/2)    R*sin(2/3*pi/2)   -W/2;
    factor*R/4+R*cos(2/3*pi/2)     R*sin(2/3*pi/2)   W/2;];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+R*cos(-1/3*pi/2)    R*sin(-1/3*pi/2)   -W/2;
    factor*R/4+R*cos(-1/3*pi/2)     R*sin(-1/3*pi/2)   W/2;
    factor*R/4+R*cos(-2/3*pi/2)    R*sin(-2/3*pi/2)   -W/2;
    factor*R/4+R*cos(-2/3*pi/2)     R*sin(-2/3*pi/2)   W/2;];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-R*cos(1/3*pi/2)    -R*sin(1/3*pi/2)   -W/2;
    factor*R/4+-R*cos(1/3*pi/2)     -R*sin(1/3*pi/2)   W/2;
    factor*R/4+-R*cos(2/3*pi/2)    -R*sin(2/3*pi/2)   -W/2;
    factor*R/4+-R*cos(2/3*pi/2)     -R*sin(2/3*pi/2)   W/2;];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-R*cos(-1/3*pi/2)    -R*sin(-1/3*pi/2)   -W/2;
    factor*R/4+-R*cos(-1/3*pi/2)     -R*sin(-1/3*pi/2)   W/2;
    factor*R/4+-R*cos(-2/3*pi/2)    -R*sin(-2/3*pi/2)   -W/2;
    factor*R/4+-R*cos(-2/3*pi/2)     -R*sin(-2/3*pi/2)   W/2;];

% second round
node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    R*cos(1/3*pi/2)   R*sin(1/3*pi/2);
    factor*R/4+W/2     R*cos(1/3*pi/2)   R*sin(1/3*pi/2);
    factor*R/4+-W/2    R*cos(2/3*pi/2)   R*sin(2/3*pi/2);
    factor*R/4+W/2     R*cos(2/3*pi/2)   R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    R*cos(-1/3*pi/2)   R*sin(-1/3*pi/2);
    factor*R/4+W/2     R*cos(-1/3*pi/2)   R*sin(-1/3*pi/2);
    factor*R/4+-W/2    R*cos(-2/3*pi/2)   R*sin(-2/3*pi/2);
    factor*R/4+W/2     R*cos(-2/3*pi/2)   R*sin(-2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    -R*cos(1/3*pi/2)   -R*sin(1/3*pi/2);
    factor*R/4+W/2     -R*cos(1/3*pi/2)   -R*sin(1/3*pi/2);
    factor*R/4+-W/2    -R*cos(2/3*pi/2)   -R*sin(2/3*pi/2);
    factor*R/4+W/2     -R*cos(2/3*pi/2)   -R*sin(2/3*pi/2);];

node.coordinates_mat=[node.coordinates_mat;
    factor*R/4+-W/2    -R*cos(-1/3*pi/2)   -R*sin(-1/3*pi/2);
    factor*R/4+W/2     -R*cos(-1/3*pi/2)   -R*sin(-1/3*pi/2);
    factor*R/4+-W/2    -R*cos(-2/3*pi/2)   -R*sin(-2/3*pi/2);
    factor*R/4+W/2     -R*cos(-2/3*pi/2)   -R*sin(-2/3*pi/2);];

% spring connector
node.coordinates_mat=[node.coordinates_mat;
   (factor-1)*R/4+R   -W/2   -W/2;
   (factor-1)*R/4+R   -W/2   W/2;
   (factor-1)*R/4+R   W/2    -W/2;
   (factor-1)*R/4+R   W/2    W/2];

tempNodeNum=size(node.coordinates_mat);
tempNodeNum=tempNodeNum(1);



%% Define assembly
assembly=Assembly_MembraneWithBar;
cst=Vec_Elements_CST;
rotSpr=Vec_Elements_RotSprings_4N;
t2t=CD_Elements_T2T_Contact;
bar=CD_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.t2t=t2t;
assembly.rotSpr=rotSpr;
assembly.bar=bar;

%% Define Plotting Functions
plots=Plot_MembraneWithBar;
plots.assembly=assembly;
plots.displayRange=[-0.2; 0.5; -0.2; 0.2; -0.2; 0.2];

plots.Plot_Shape_NodeNumber;


%% Define Triangle
% First ball

cst.node_ijk_mat=[cst.node_ijk_mat;
    1 2 3;
    2 3 4;
    2 4 27;
    4 27 28;
    27 28 26;
    25 26 27;
    22 26 25;
    24 22 26;
    1 3 40;
    1 39 40;
    37 38 39;
    38 39 40;
    18 20 37;
    20 37 38;
    17 18 20;
    17 19 20;
    17 34 33;
    34 19 17;
    33 34 35;
    35 36 34;
    5 7 35;
    35 36 7;
    5 6 7;
    6 7 8;
    6 8 31;
    31 32 8;
    31 32 29;
    30 29 32;
    29 30 21;
    21 23 30;
    21 24 23;
    21 22 24;
    9 10 11;
    11 10 12;
    12 10 53;
    53 54 12;
    53 54 56;
    56 55 53;
    56 55 2;
    1 2 55;
    5 6 51;
    51 52 6;
    9 11 49;
    49 50 11;
    49 50 51;
    51 52 50;
    3 4 43;
    43 44 4;
    41 42 43;
    43 44 42;
    14 16 41;
    41 42 16;
    13 14 15;
    14 16 15;
    13 15 46;
    45 46 13;
    45 46 47;
    47 48 46;
    7 8 47;
    47 48 8;
    ];


%% Second ball

cst.node_ijk_mat=[cst.node_ijk_mat;
    66 67 68;
    65 66 67;
    65 66 90;
    89 90 65;
    89 90 91;
    91 92 90;
    91 92 74;
    73 74 91;
    73 74 75;
    74 75 76;
    75 76 87;
    76 87 88;
    87 88 86;
    87 85 86;
    85 86 70;
    69 70 85;
    69 70 71;
    71 72 70;
    71 72 82;
    81 82 71;
    81 83 82;
    82 83 84;
    83 84 79;
    79 80 84;
    77 78 79;
    78 79 80;
    77 78 95;
    95 96 78;
    93 94 96;
    95 96 93;
    67 93 94;
    67 68 94;
    57 58 59;
    58 59 60;
    57 59 112;
    111 112 57;
    109 110 111;
    111 112 110;
    74 76 109;
    109 110 76;
    73 75 105;
    105 106 75;
    105 106 107;
    107 108 106;
    107 108 61;
    61 63 108;
    61 62 63;
    62 63 64;
    62 64 103;
    103 104 64;
    103 104 101;
    101 102 104;
    77 79 101;
    101 102 79;
    78 80 98;
    97 98 78;
    97 98 99;
    99 100 98;
    60 100 99;
    58 60 99;];

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);

% material properties
t=0.0004;
E=2.5*10^9;
cst.t_vec=t*ones(cstNum,1); % Thickness 
cst.E_vec=E*ones(cstNum,1); % Young's Modulus
cst.v_vec=0.2*ones(cstNum,1); % Poisson's Ratio

plots.Plot_Shape_CSTNumber;
plots.Plot_Shape_NodeNumber;

%% Define Rotational Spring
% First Ball
rotSpr.node_ijkl_mat=[rotSpr.node_ijkl_mat;
    1 2 3 4;
    3 2 4 27;
    2 4 27 28;
    4 27 28 26;
    28 27 26 25
    27 26 25 22;
    25 26 22 24;
    26 22 24 21;
    22 24 21 23;
    24 21 23 30;
    23 21 30 29;
    21 30 29 32;
    30 29 32 31;
    29 32 31 8;
    32 31 8 6;
    31 8 6 7;
    8 6 7 5;
    6 7 5 35;
    5 7 35 36;
    7 35 36 34;
    36 35 34 33;
    35 33 34 17;
    33 34 17 19;
    34 17 19 20;
    19 17 20 18;
    17 18 20 37;
    18 20 37 38;
    20 37 38 39;
    37 38 39 40;
    38 39 40 1;
    39 40 1 3;
    40 1 3 2;
    2 4 3 43;
    3 4 43 44;
    4 43 44 42;
    44 43 42 41;
    43 42 41 16;
    42 41 16 14;
    41 16 14 15;
    16 14 15 13;
    14 15 13 46;
    15 13 46 45;
    13 46 45 47;
    45 46 47 48;
    46 47 48 8;
    48 47 8 7;
    47 8 7 6;
    7 6 5 51;
    5 6 51 52;
    8 51 52 50;
    52 51 50 49;
    51 50 49 11;
    50 49 11 9;
    49 11 9 10;
    9 11 10 12;
    11 10 12 53;
    10 12 53 54;
    12 53 54 56;
    54 53 56 55;
    53 56 55 2;
    56 55 2 1;
    55 2 1 3;
    ];

%% Second Ball
rotSpr.node_ijkl_mat=[
    rotSpr.node_ijkl_mat;
    57 59 58 60;
    59 58 60 99;
    58 60 99 100;
    60 99 100 98;
    100 99 98 97;
    99 98 97 78;
    97 98 78 80;
    98 78 80 79;
    80 78 79 77;
    78 79 77 101;
    77 79 101 102;
    79 101 102 104;
    102 101 104 103;
    101 104 103 64;
    104 103 64 62;
    103 64 62 63;
    64 62 63 61;
    62 63 61 106;
    63 61 108 107;
    61 108 107 106;
    108 107 106 105;
    107 106 105 75;
    106 105 75 73;
    105 75 73 74;
    73 75 74 76;
    75 74 76 109;
    74 76 109 110;
    76 109 110 111;
    109 110 111 112;
    110 111 112 57;
    111 112 57 59;
    112 57 59 58;
    74 75 76 87;
    75 76 87 88;
    76 87 88 86;
    88 87 86 85;
    87 86 85 70;
    86 85 70 69;
    85 70 69 71;
    69 70 71 72;
    70 71 72 82;
    72 71 82 81;
    71 82 81 83;
    81 82 83 84;
    82 83 84 79;
    83 84 79 80;
    84 79 80 78;
    79 77 78 95;
    77 78 95 96;
    78 95 96 93;
    95 96 93 94;
    96 93 94 67;
    93 94 67 68;
    94 67 68 66;
    68 67 66 65;
    67 66 65 90;
    66 65 90 89;
    65 90 89 91;
    89 90 91 92;
    90 91 92 74;
    92 91 74 73;
    91 74 73 75;
    ];

rotNum=size(rotSpr.node_ijkl_mat);
rotNum=rotNum(1);

% Find the rotational spring stiffness
L=2*pi*R/16;
rotk=E*t^3*W/12/L;

rotSpr.rot_spr_K_vec=rotk*ones(rotNum,1);
rotSpr.theta_stress_free_vec=pi*ones(rotNum,1);
plots.Plot_Shape_SprNumber;


%% Define Triangle to Triangle Penertration Prevention
t2t.d0=0.01;
t2t.delta=10^-6;
t2t.k_contact=10;

t2t.tri_ijk_mat=cst.node_ijk_mat;
t2t.group_number=[t2t.group_number;
              ones(cstNum/2,1);];

t2t.group_number=[t2t.group_number;
              ones(cstNum/2,1);];

% Only one contact group will 
% deactivate the penetration prevention





%% Set up bar
bar.node_ij_mat=[69 113;
                 70 114;
                 71 115;
                 72 116];
bar.E_vec=ones(4,1);
bar.A_vec=0.02*ones(4,1);


% set up ploting color
plots.colorNum=[ones(cstNum/2,1);
              2*ones(cstNum/2,1);];

assembly.Initialize_Assembly;



%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nr.supp=[9 1 1 1;
     10 1 1 1;
     11 1 1 1;
     12 1 1 1;
     69 0 1 1;
     70 0 1 1;
     71 0 1 1;
     72 0 1 1;
     113 1 1 1;
     114 1 1 1;
     115 1 1 1;
     116 1 1 1;];

force=0.004;
step=20;

Uhis=[];
for k=1:step
    nr.load=[69 force*k 0 0;
             70 force*k 0 0;
             71 force*k 0 0;
             72 force*k 0 0];
    
    nr.increStep=1;
    nr.iterMax=5;
    nr.tol=1*10^-2;

    Uhis(k,:,:)=squeeze(nr.Solve());
end

toc
plots.Plot_DeformedShape(squeeze(Uhis(end,:,:)))
plots.Plot_DeformedShape(zeros(size(squeeze(Uhis(end,:,:)))))

plots.fileName='InterlockedBalls_NoContact.gif';
plots.Plot_DeformedHis(Uhis)

RefUHis=squeeze(Uhis(:,[69 70 71 72],1));
RefUHis=mean(RefUHis,2);

forceHis=(1:step)*force*4;
forceHis=forceHis';

figure
plot([0,RefUHis'],[0,forceHis']);
xlabel('displacement (m)')
ylabel('force (N)')


