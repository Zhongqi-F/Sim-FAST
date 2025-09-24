clear all;
clc;
close all;

%% Initialize the scissor 
N=8;
H=1; % (m)
L=1; % (m)

barA=0.01; 
barE=2*10^9; % (Pa)
panel_E=200*10^9;
panel_t=0.05;
panel_v=0.2;

I=1/12*0.01^4; 
barL=sqrt(H^2+L^2);
kspr=3*barE*I/barL;
node=Elements_Nodes();

%% Define the nodal coordinates
for i=1:N
    node.coordinates_mat=[node.coordinates_mat;
        L*(i-1), 0, 0;
        L*(i-1), L, 0;
        L*(i-1), 0, L;
        L*(i-1), L, L;

        0.25*L+L*(i-1), 0, 0.25*L;
        0.25*L+L*(i-1), L, 0.25*L;
        0.25*L+L*(i-1), 0, 0.75*L;
        0.25*L+L*(i-1), L, 0.75*L; %8

        0.75*L+L*(i-1), 0, 0.25*L;
        0.75*L+L*(i-1), L, 0.25*L;
        0.75*L+L*(i-1), 0, 0.75*L;
        0.75*L+L*(i-1), L, 0.75*L;

        (i-1)*L+L/2, 0, L/2;
        (i-1)*L+L/2, L, L/2;
        (i-1)*L+L/2, 0, 0;
        (i-1)*L+L/2, L, 0; %16
           
        (i-1)*L+L/2, 0, L;
        (i-1)*L+L/2, L, L;
        (i-1)*L+L/2, L/2, L;
        ];
end

node.coordinates_mat=[node.coordinates_mat;
        L*N, 0, 0;
        L*N, L, 0;
        L*N, 0, L;
        L*N, L, L;
        ];

%% Define assembly
assembly=Assembly_Scissor_Bridge(); 
cst=Vec_Elements_CST;
rotSpr3N=CD_Elements_RotSprings_3N; 
rotSpr4N=Vec_Elements_RotSprings_4N;
bar=CD_Elements_Bars;

assembly.cst=cst;
assembly.node=node;
assembly.bar=bar;
assembly.rot_spr_3N=rotSpr3N;
assembly.rot_spr_4N=rotSpr4N;

%% Define Plotting Functions
plots=Plot_Scissor_Bridge; 
plots.assembly=assembly;
plots.displayRange=[-0.5; N+0.5; -0.5; 1.5; -0.5; 1.5]; 
plots.Plot_Shape_Node_Number;

%% Define Triangle
for i=1:N
cst.node_ijk_mat=[cst.node_ijk_mat;
    19*(i-1)+1  19*(i-1)+2  19*(i-1)+15;
    19*(i-1)+2  19*(i-1)+15 19*(i-1)+16;
    19*(i-1)+16 19*(i-1)+15 19*(i-1)+20;
    19*(i-1)+20 19*(i-1)+16 19*(i-1)+21;
];
end

cstNum=size(cst.node_ijk_mat);
cstNum=cstNum(1);
cst.t_vec=panel_t*ones(cstNum,1); % m
cst.E_vec=panel_E*ones(cstNum,1); % Pa
cst.v_vec=panel_v*ones(cstNum,1); %

plots.Plot_Shape_CST_Number;

%% Define bar
for i=1:N
    bar.node_ij_mat=[
        bar.node_ij_mat;
        19*(i-1)+1   19*(i-1)+5;
        19*(i-1)+5   19*(i-1)+13;
        19*(i-1)+13  19*(i-1)+11;
        19*(i-1)+11  19*(i-1)+22;

        19*(i-1)+3   19*(i-1)+7;
        19*(i-1)+7   19*(i-1)+13;
        19*(i-1)+13  19*(i-1)+9;
        19*(i-1)+9   19*(i-1)+20;

        19*(i-1)+2   19*(i-1)+6;
        19*(i-1)+6   19*(i-1)+14;
        19*(i-1)+14  19*(i-1)+12;
        19*(i-1)+12  19*(i-1)+23;

        19*(i-1)+4   19*(i-1)+8;
        19*(i-1)+8   19*(i-1)+14;
        19*(i-1)+14  19*(i-1)+10;
        19*(i-1)+10  19*(i-1)+21;
        
        % Start of top
        19*(i-1)+3   19*(i-1)+4;
        19*(i-1)+3   19*(i-1)+17;
        19*(i-1)+3   19*(i-1)+19;        
        19*(i-1)+4   19*(i-1)+19;
        19*(i-1)+4   19*(i-1)+18;

        19*(i-1)+17  19*(i-1)+19;
        19*(i-1)+18  19*(i-1)+19;
        
        19*(i-1)+18  19*(i-1)+23;               
        19*(i-1)+19   19*(i-1)+23;
        19*(i-1)+17  19*(i-1)+22; 
        19*(i-1)+19  19*(i-1)+22;

        % add new bars
        19*(i-1)+1    19*(i-1)+2;
        19*(i-1)+15   19*(i-1)+16;
        19*(i-1)+1    19*(i-1)+15;
        19*(i-1)+15   19*(i-1)+20;
        19*(i-1)+2    19*(i-1)+16;
        19*(i-1)+16   19*(i-1)+21;
        19*(i-1)+20   19*(i-1)+21;
        19*(i-1)+2    19*(i-1)+15;
        19*(i-1)+16   19*(i-1)+20;
        ];
end

bar.node_ij_mat=[
        bar.node_ij_mat;
        19*N+3    19*N+4;
        ];

barNum=size(bar.node_ij_mat);
barNum=barNum(1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);
plots.Plot_Shape_Bar_Number();

%% Define 3 Node Rotational Spring
for i=1:N    
    rotSpr3N.node_ijk_mat = [rotSpr3N.node_ijk_mat;
         19*(i-1)+1  19*(i-1)+5  19*(i-1)+13;
         19*(i-1)+5  19*(i-1)+13 19*(i-1)+11;
         19*(i-1)+13 19*(i-1)+11 19*(i-1)+22;

         19*(i-1)+3  19*(i-1)+7  19*(i-1)+13;
         19*(i-1)+7  19*(i-1)+13 19*(i-1)+9;
         19*(i-1)+13 19*(i-1)+9  19*(i-1)+20;

         19*(i-1)+2  19*(i-1)+6  19*(i-1)+14;
         19*(i-1)+6  19*(i-1)+14 19*(i-1)+12;
         19*(i-1)+14 19*(i-1)+12 19*(i-1)+23;

         19*(i-1)+4  19*(i-1)+8  19*(i-1)+14;
         19*(i-1)+8  19*(i-1)+14 19*(i-1)+10;
         19*(i-1)+14 19*(i-1)+10 19*(i-1)+21;
         ];
end

rotNum=size(rotSpr3N.node_ijk_mat,1);
rotSpr3N.rot_spr_K_vec=kspr*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;



%% Set up four node rotational spring
for i=1:N    
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
         19*(i-1)+1   19*(i-1)+2   19*(i-1)+15  19*(i-1)+16;
         19*(i-1)+2   19*(i-1)+15  19*(i-1)+16  19*(i-1)+20;
         19*(i-1)+15  19*(i-1)+16  19*(i-1)+20  19*(i-1)+21;

         19*(i-1)+3  19*(i-1)+17  19*(i-1)+19  19*(i-1)+22;
         19*(i-1)+4  19*(i-1)+18  19*(i-1)+19  19*(i-1)+23;

         19*(i-1)+3  19*(i-1)+4  19*(i-1)+19  19*(i-1)+18;
         19*(i-1)+4  19*(i-1)+3  19*(i-1)+19  19*(i-1)+17;
         19*(i-1)+18  19*(i-1)+19  19*(i-1)+23  19*(i-1)+22;
         19*(i-1)+17  19*(i-1)+22  19*(i-1)+19  19*(i-1)+23;
         ];
end

rotNum4N=size(rotSpr4N.node_ijkl_mat,1);
rotSpr4N.rot_spr_K_vec=10000*ones(rotNum4N,1);
plots.Plot_Shape_RotSpr_4N_Number;
assembly.Initialize_Assembly;


%% Calculate Self-weight of the Bridge
rho_steel=7850;         % Density of steel in kg/m^3
g=9.81;                 % Gravitational acceleration in m/s^2

% Bar elements
A_bar=barA;             % Cross-sectional area of bars in m^2 (matches bar.A_vec)

% Calculate total length of all bars
L_total=0;
barNodeMat=bar.node_ij_mat;
coords=node.coordinates_mat;

for i=1:size(barNodeMat,1)
    n1=barNodeMat(i,1);
    n2=barNodeMat(i,2);
    p1=coords(n1,:);
    p2=coords(n2,:);
    len=norm(p1-p2);
    L_total=L_total+len;
end

W_bar=A_bar*L_total*rho_steel*g;   % Total bar weight (Newtons)

% CST panel elements
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

t_cst=panel_t; 
W_cst=A_cst_total*t_cst*rho_steel*g; % Total CST panel weight (Newtons)

% Total self-weight
W_total=W_bar+W_cst;

% Output results
fprintf('-----------------------------\n');
fprintf('Total length of all bars: %.3f m\n', L_total);
fprintf('Total area of all CST panels: %.3f m^2\n', A_cst_total);
fprintf('Total bar weight: %.2f N\n', W_bar);
fprintf('Total CST panel weight: %.2f N\n', W_cst);
fprintf('Total self-weight of the bridge: %.2f N\n', W_total);
fprintf('-----------------------------\n');


%% Set up solver   
nr = Solver_NR_Loading;
nr.assembly = assembly;

nr.supp = [1    1 1 1;
           2    1 1 1;
           38*N/2+1    1 1 1; 
           38*N/2+2    1 1 1; 
           ];


force=1000;
step=10;

nr.increStep = 10;
nr.load=[38*N/4+1 0 0 -force/2; 
         38*N/4+2 0 0 -force/2]; 

nr.iterMax = 20;
nr.tol = 10^-5;

Uhis = nr.Solve();  


U_end = squeeze(Uhis(end, :, :));  
plots.Plot_Deformed_Shape(20*U_end);



% Failure load computation
truss_strain=bar.Solve_Strain(node,squeeze(Uhis(end,:,:)));
internal_force=(truss_strain).*(bar.E_vec).*(bar.A_vec);

% Find the maximum bar internal force
[maxBarForce,index]=max(abs(internal_force));
% [maxBarForce,index]=max(internal_force);


% Find failure force for the bar
sigma_u=300*10^6;
barFailureForce=sigma_u*barA;


% Find total bar length
barLtotal=sum(bar.L0_vec);


% Find Stiffness
Uaverage=-mean(squeeze(Uhis(end,[38*N/4+1,38*N/4+2],3)));
Kstiff=step*force/Uaverage;


% Plot failure stress
bar_stress=(truss_strain).*(bar.E_vec)*barFailureForce/maxBarForce;
plots.Plot_Shape_Bar_Stress(bar_stress)


% Find the relationship betweeen the bar internal forces and load
loadatfail=force*step*barFailureForce/maxBarForce;
fprintf('Failure load is %d kN \n',  loadatfail/1000);
fprintf('Total bar length is %d m \n',  barLtotal);
fprintf('Stiffness is %d N/m \n',  Kstiff);
fprintf('Total bars: %d\n', barNum);