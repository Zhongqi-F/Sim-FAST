clear all;
clc;
close all;
tic

%% Define the Rolling Bridge Geometry
% Height of the bridge
H=2; % meter

% Width of the bridge
W=2; % meter

% Length of the section
L=2; % meter

% Number of Sections
N=8;

% Load the deformation history
UhisNew=load('RollingUhis.mat');
UhisNew=UhisNew.Uhis;
DepRate=1; % 1 is fully deployed, 0 is compact
DepStep = max(1, int32((1 - DepRate) * 800));

% The cross section of this bridge is:
% HSS 4X3X5/16 A500 Grade C Fy=50ksi/200 GPa
barA=0.0023;
barE=2*10^11;

% We will use the weak axis
I=1.88*10^-6;

% We assume a soft panel so that only the truss is taking global load
% Thus, panel Young's modulus is 200 MPa
panel_E=2*10^8;
panel_t=0.01;
panel_v=0.3;

% We assume that the actuator bar has similar material property
% and stiffness as normal truss
activeBarE=2*10^11;

% Initialize Elements and Assembly
node=Elements_Nodes;
bar=Vec_Elements_Bars;
actBar=CD_Elements_Bars;
cst=Vec_Elements_CST;
rot_spr_4N=Vec_Elements_RotSprings_4N;

assembly=Assembly_Rolling_Bridge;
assembly.node=node;
assembly.bar=bar;
assembly.actBar=actBar;
assembly.cst=cst;
assembly.rot_spr_4N=rot_spr_4N;


%% Define the nodal coordinates
node.coordinates_mat=[];
node.coordinates_mat=[
    node.coordinates_mat;
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


%% Define Plotting Functions
plots=Plot_Rolling_Bridge();
plots.assembly=assembly;

plots.displayRange=[-2;18;-1;3;-1;14];
plots.viewAngle1=20;
plots.viewAngle2=20;

plots.Plot_Shape_Node_Number()


%% Define triangle elements
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

% Add top bars for stability
bar.node_ij_mat=[bar.node_ij_mat;
    2 6;
    7 8;
    ];

for i=1:N-1
    bar.node_ij_mat=[bar.node_ij_mat;
        10+(i-2)*6,  12+(i-2)*6;
        13+(i-2)*6,  14+(i-2)*6;
        16+(i-2)*6,  18+(i-2)*6;
        ];
end

barNum=size(bar.node_ij_mat,1);
bar.A_vec=barA*ones(barNum,1);
bar.E_vec=barE*ones(barNum,1);

plots.Plot_Shape_Bar_Number();


%% Define how actuator bars are connected
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


%% Define four node rotational springs

% Major Rot Spr
rot_spr_4N.node_ijkl_mat=[
    4 1 3 2;
    3 4 5 6;];

for i=1:N-1
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+5, 6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+10;
        6*(i-1)+9, 6*(i-1)+5, 6*(i-1)+11, 6*(i-1)+12;];
end

% Floor Rot Spr
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    1 3 4 5];

for i=1:N-1
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+5, 6*(i-1)+11;];
end

% Secondary Rot Spring
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    1 2 3 7;
    4 5 6 8;];

for i=1:N-2
    rot_spr_4N.node_ijkl_mat=[
        rot_spr_4N.node_ijkl_mat;
        6*(i-1)+7, 6*(i-1)+3, 6*(i-1)+10, 6*(i-1)+9;
        6*(i-1)+3, 6*(i-1)+9, 6*(i-1)+10, 6*(i-1)+13;
        6*(i-1)+8, 6*(i-1)+5, 6*(i-1)+12, 6*(i-1)+11;
        6*(i-1)+5, 6*(i-1)+11, 6*(i-1)+12, 6*(i-1)+14;];
end

i=N-1;
rot_spr_4N.node_ijkl_mat=[
    rot_spr_4N.node_ijkl_mat;
    6*(i-1)+7, 6*(i-1)+3, 6*(i-1)+10, 6*(i-1)+9;
    6*(i-1)+8, 6*(i-1)+5, 6*(i-1)+12, 6*(i-1)+11;
    ];

rotNum=size(rot_spr_4N.node_ijkl_mat,1);
rot_spr_4N.rot_spr_K_vec=(10^6)*ones(rotNum,1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;


%% Initialize the entire assembly
node.coordinates_mat=node.coordinates_mat+squeeze(UhisNew(DepStep,:,:));
assembly.Initialize_Assembly();


%% Calculate Self-weight of the Bridge

% Steel properties
rho_steel=7850;         % Density of steel (kg/m^3)
g=9.81;                 % Gravitational acceleration (m/s^2)

% Bar cross-section area
A_bar=barA;             % m^2, should match bar.A_vec

% Total length of all normal bar elements
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

% Total length of all actuator bar elements
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

% Total bar length (all bars)
L_total=L_bar_total+L_actbar_total;

% Total weight of all bars (Newtons)
W_bar=A_bar*L_total*rho_steel*g;


%% Set up solver
nr=Solver_NR_Loading;
nr.assembly=assembly;

nodeNum=size(node.coordinates_mat,1);
nodeNumVec=(1:nodeNum)';

% Set up the support of this bridge
nr.assembly=assembly;
nr.supp=[nodeNumVec,zeros(nodeNum,1),zeros(nodeNum,1),zeros(nodeNum,1)];
nr.supp(1,2:4)=ones(1,3);
nr.supp(4,2:4)=ones(1,3);
nr.supp(3,2:4)=ones(1,3);
nr.supp(5,2:4)=ones(1,3);


% -----------------------------------------------------------------------
%  AASHTO LRFD Material & Section Properties
%  (defined once outside the loop — these never change)
% -----------------------------------------------------------------------
Fy  = 345e6;    % Yield strength (Pa), ASTM A500 Gr.C
Fu  = 427e6;    % Tensile strength (Pa), A500 Gr.C
E   = barE;     % Elastic modulus (Pa)

% Resistance factors — AASHTO LRFD Art. 6.5.4.2
phi_ty = 0.95;  % Tension yielding (gross section)
phi_uf = 0.80;  % Tension fracture (net section)
phi_c  = 0.95;  % Axial compression

% Net section parameters (welded connection)
% AASHTO LRFD Art. 6.8.2.1
Ag = barA;
An = barA;   % No bolt holes -> net area = gross area
U  = 1.0;    % Shear lag factor
Rp = 1.0;    % Hole reduction factor (drilled/reamed holes)

% Radius of gyration (weak axis governs)
r_val = sqrt(I / barA);

% Effective length factor
K = 1.0;

% -----------------------------------------------------------------------
%  LOCAL BUCKLING CHECK — AASHTO LRFD Art. 6.9.4.2
%  Section: HSS 4x3x5/16, cold-formed, A500 Gr.C
%  Computed once — section geometry does not change with load steps
% -----------------------------------------------------------------------
bt = 7.31;   % b/t ratio of wider face  (flange), from AISC tables
ht = 10.7;   % h/t ratio of narrower face (web),  from AISC tables

% Limiting slenderness for uniformly compressed plate elements in HSS
% AASHTO LRFD Art. 6.9.4.2.1, Table 6.9.4.2.1-1
lambda_r = 1.92 * sqrt(E / Fy);

local_buckle_bt_pass = (bt <= lambda_r);
local_buckle_ht_pass = (ht <= lambda_r);
local_buckle_pass    = local_buckle_bt_pass && local_buckle_ht_pass;

fprintf('--- Local Buckling Check (AASHTO LRFD Art. 6.9.4.2) ---\n');
fprintf('  Limiting slenderness lambda_r = 1.92*sqrt(E/Fy) = %.2f\n', lambda_r);
fprintf('  b/t = %.2f  ->  %s\n', bt, sel_str(local_buckle_bt_pass, 'OK', 'FAIL'));
fprintf('  h/t = %.2f  ->  %s\n', ht, sel_str(local_buckle_ht_pass, 'OK', 'FAIL'));
if local_buckle_pass
    fprintf('  Section is non-slender (local buckling OK)\n');
else
    fprintf('  WARNING: Section FAILS local buckling slenderness limit\n');
end
fprintf('--------------------------------------------------------\n');


maxDCR_val  = NaN;
maxDCR_idx  = NaN;
DCR         = NaN;
phi_Rn      = NaN;
modeStr     = {};
slender_chk = {};
passYN      = [];

% -----------------------------------------------------------------------
%  LOAD STEP LOOP
% -----------------------------------------------------------------------
for i = 1:50

    % Nonlinear solver settings
    nr.increStep = 1;
    nr.iterMax   = 50;
    nr.tol       = 1e-5;

    % Apply self weight incrementally
    nodeNum = size(node.coordinates_mat, 1);
    force   = W_bar / nodeNum / 50 * i;

    nr.load = [(1:nodeNum)'  zeros(nodeNum,1) ...
                zeros(nodeNum,1)  -force*ones(nodeNum,1)];

    % Solve
    Uhis = nr.Solve;

    % -------------------------------------------------------------------
    %  Evaluate if member is failing (AASHTO LRFD)
    % -------------------------------------------------------------------
    U_end = squeeze(Uhis(end,:,:));   % [nodeNum x 3]

    % Strain and internal force in each bar
    truss_strain   = bar.Solve_Strain(node, U_end);
    internal_force = truss_strain .* (bar.E_vec) .* (bar.A_vec);

    % Replace NaN/Inf from degenerate (near-zero-length) members with 0
    internal_force(~isfinite(internal_force)) = 0;
    truss_strain(~isfinite(truss_strain))     = 0;

    barNum  = numel(internal_force);
    L0_vec  = bar.L0_vec(:);
    Lc      = K .* L0_vec;

    % Output arrays (reset each load step)
    passYN      = false(barNum, 1);
    DCR         = NaN(barNum, 1);
    phi_Rn      = NaN(barNum, 1);
    modeStr     = cell(barNum, 1);
    slender_chk = cell(barNum, 1);

    % Member-by-member check
    for k = 1:barNum

        Pu  = internal_force(k);
        Lck = Lc(k);

        % --- Guard against degenerate members (zero-length or NaN force) --
        if isnan(Pu) || isnan(Lck) || Lck < 1e-9
            DCR(k)         = NaN;
            phi_Rn(k)      = NaN;
            modeStr{k}     = 'SKIP (degenerate member)';
            slender_chk{k} = 'SKIP';
            passYN(k)      = true;  % exclude from governing check
            continue
        end

        % --- Global slenderness check ------------------------------------
        slender_ratio = Lck / r_val;

        if Pu < 0
            % Compression: KL/r <= 120  (AASHTO LRFD Art. 6.9.3)
            slender_limit = 120;
            if slender_ratio > slender_limit
                slender_chk{k} = sprintf('FAIL KL/r=%.1f > %d', slender_ratio, slender_limit);
            else
                slender_chk{k} = sprintf('OK   KL/r=%.1f <= %d', slender_ratio, slender_limit);
            end
        else
            % Tension: L/r <= 200  (AASHTO LRFD Art. 6.8.4)
            slender_limit = 200;
            if slender_ratio > slender_limit
                slender_chk{k} = sprintf('FAIL L/r=%.1f > %d', slender_ratio, slender_limit);
            else
                slender_chk{k} = sprintf('OK   L/r=%.1f <= %d', slender_ratio, slender_limit);
            end
        end

        % --- Resistance calculation --------------------------------------
        if Pu >= 0
            % TENSION MEMBER (AASHTO LRFD Art. 6.8.2.1)
            phiRn_ty  = phi_ty * Fy * Ag;           % Yielding of gross section
            phiRn_uf  = phi_uf * Fu * An * Rp * U;  % Fracture of net section
            phi_Rn(k) = min(phiRn_ty, phiRn_uf);

            if phi_Rn(k) == phiRn_ty
                modeStr{k} = 'Tension Yielding (Art.6.8.2.1)';
            else
                modeStr{k} = 'Tension Fracture (Art.6.8.2.1)';
            end

            DCR(k) = abs(Pu) / phi_Rn(k);

        else
            % COMPRESSION MEMBER (AASHTO LRFD Art. 6.9.4.1)
            Fe = (pi^2 * E) / (slender_ratio^2);    % Elastic buckling stress

            if (Fy / Fe) <= 2.25
                Fcr        = 0.658^(Fy / Fe) * Fy;
                modeStr{k} = 'Compression Inelastic Buckling (Art.6.9.4.1)';
            else
                Fcr        = 0.877 * Fe;
                modeStr{k} = 'Compression Elastic Buckling (Art.6.9.4.1)';
            end

            phi_Rn(k) = phi_c * Fcr * Ag;
            DCR(k)    = abs(Pu) / phi_Rn(k);
        end

        % --- Pass/Fail ---------------------------------------------------
        slender_pass = ~contains(slender_chk{k}, 'FAIL');

        if Pu < 0
            % Compression: DCR + global slenderness + local buckling
            passYN(k) = (DCR(k) <= 1.0) && slender_pass && local_buckle_pass;
        else
            % Tension: DCR + global slenderness only
            passYN(k) = (DCR(k) <= 1.0) && slender_pass;
        end

    end  % end member loop

    % -------------------------------------------------------------------
    %  Safety check — break loop if any member fails
    % -------------------------------------------------------------------
    [maxDCR_val, maxDCR_idx] = max(DCR, [], 'omitnan');

    if all(passYN)
        fprintf('Step %2d/50 : All Truss Members Safe (AASHTO LRFD)\n', i);
    else
        fprintf('Step %2d/50 : Member Failure Detected (AASHTO LRFD)\n', i);
        break
    end

end  % end load step loop


% -----------------------------------------------------------------------
%  Linear extrapolation capacity estimate
% -----------------------------------------------------------------------
SF_linear      = 1.0 / maxDCR_val;
W_capacity_lin = W_bar * SF_linear;

fprintf('\n========================================\n');
fprintf('  BRIDGE CAPACITY ANALYSIS (AASHTO LRFD)\n');
fprintf('========================================\n');
fprintf('  Self-weight W_bar        = %.2f N\n',   W_bar);
fprintf('  Max DCR at W_bar         = %.4f\n',     maxDCR_val);
fprintf('  Safety Factor            = %.4f\n',     SF_linear);
fprintf('  Bridge Capacity (approx) = %.2f N\n',   W_capacity_lin);
fprintf('  (= %.4f x self-weight)\n',              SF_linear);
fprintf('========================================\n');


% -----------------------------------------------------------------------
%  Find average tip deflection (from last solved step)
% -----------------------------------------------------------------------
Uaverage=-mean(squeeze(Uhis(end,[45,47],3)));


% -----------------------------------------------------------------------
%  Output results (AASHTO LRFD)
% -----------------------------------------------------------------------
fprintf('-----------------------------\n');
fprintf('Total length of all bars                  : %.2f m\n',  L_total);
fprintf('Total bar weight                          : %.2f N\n',  W_bar);
fprintf('Local buckling b/t = %.2f                 : %s\n',      bt, sel_str(local_buckle_bt_pass, 'OK', 'FAIL'));
fprintf('Local buckling h/t = %.2f                 : %s\n',      ht, sel_str(local_buckle_ht_pass, 'OK', 'FAIL'));
fprintf('Maximum DCR (Demand/Capacity Ratio)       : %.4f\n',    maxDCR_val);
fprintf('Governing bar index                       : %d\n',      maxDCR_idx);
fprintf('Governing limit state                     : %s\n',      modeStr{maxDCR_idx});
fprintf('Governing bar slenderness                 : %s\n',      slender_chk{maxDCR_idx});
fprintf('Design resistance phiRn (governing bar)   : %.2f N\n',  phi_Rn(maxDCR_idx));
fprintf('Tip deflection                            : %.6f m\n',  Uaverage);
fprintf('-----------------------------\n');

% Plot bar stress
truss_stress = truss_strain .* (bar.E_vec);
plots.Plot_Shape_Bar_Stress(truss_stress);

% Plot failed bars
plots.Plot_Shape_Bar_Failure(passYN);

% Plot deformed shape
plots.Plot_Deformed_Shape(squeeze(Uhis(end,:,:)));


% -----------------------------------------------------------------------
%  Critical Member Detailed Report (AASHTO LRFD)
% -----------------------------------------------------------------------
fprintf('\n========================================\n');
fprintf('  CRITICAL MEMBER REPORT (AASHTO LRFD)  \n');
fprintf('========================================\n');

% Node indices of the critical member
n1_crit = bar.node_ij_mat(maxDCR_idx, 1);
n2_crit = bar.node_ij_mat(maxDCR_idx, 2);

% Internal force at the last solved load step
Pu_crit = internal_force(maxDCR_idx);

% Effective length and slenderness ratio
Lc_crit = K * bar.L0_vec(maxDCR_idx);
sr_crit  = Lc_crit / r_val;

fprintf('  Member index               : %d\n',       maxDCR_idx);
fprintf('  Connected nodes            : %d -- %d\n', n1_crit, n2_crit);
fprintf('  Internal force Pu          : %.2f N  (%s)\n', ...
        abs(Pu_crit), sel_str(Pu_crit >= 0, 'Tension', 'Compression'));
fprintf('  Design resistance phi*Rn   : %.2f N\n',   phi_Rn(maxDCR_idx));
fprintf('  DCR (Demand/Capacity)      : %.4f  ->  %s\n', ...
        maxDCR_val, sel_str(maxDCR_val <= 1.0, 'PASS', 'FAIL'));
fprintf('  Governing limit state      : %s\n',       modeStr{maxDCR_idx});
fprintf('  Slenderness check          : %s\n',       slender_chk{maxDCR_idx});
fprintf('  Effective length Lc        : %.4f m\n',   Lc_crit);
fprintf('  Slenderness ratio          : %.2f\n',     sr_crit);
fprintf('  Local buckling check       : %s\n',       ...
        sel_str(local_buckle_pass, 'OK', 'FAIL'));
fprintf('  Overall member status      : %s\n',       ...
        sel_str(passYN(maxDCR_idx), 'PASS', 'FAIL'));
fprintf('========================================\n');


% -----------------------------------------------------------------------
%  Visualize deployed shape at current DepRate and highlight critical bar
% -----------------------------------------------------------------------
figure;
hold on;
axis equal;
view(plots.viewAngle1, plots.viewAngle2);

coords_plot = node.coordinates_mat;
barMat      = bar.node_ij_mat;
cstMat      = cst.node_ijk_mat;

% -----------------------------
% Plot CST elements as yellow triangles
% -----------------------------
for e = 1:size(cstMat,1)
    n1 = cstMat(e,1);
    n2 = cstMat(e,2);
    n3 = cstMat(e,3);

    X = [coords_plot(n1,1), coords_plot(n2,1), coords_plot(n3,1)];
    Y = [coords_plot(n1,2), coords_plot(n2,2), coords_plot(n3,2)];
    Z = [coords_plot(n1,3), coords_plot(n2,3), coords_plot(n3,3)];

    patch(X, Y, Z, 'y', ...
        'FaceAlpha', 0.8, ...
        'EdgeColor', 'none');
end

% -----------------------------
% Plot all bars in black
% -----------------------------
for e = 1:size(barMat,1)
    n1 = barMat(e,1);
    n2 = barMat(e,2);

    x = [coords_plot(n1,1), coords_plot(n2,1)];
    y = [coords_plot(n1,2), coords_plot(n2,2)];
    z = [coords_plot(n1,3), coords_plot(n2,3)];

    plot3(x, y, z, 'k-', 'LineWidth', 1.2);
end

% -----------------------------
% Highlight critical bar in red
% -----------------------------
n1c = barMat(maxDCR_idx,1);
n2c = barMat(maxDCR_idx,2);

plot3([coords_plot(n1c,1), coords_plot(n2c,1)], ...
      [coords_plot(n1c,2), coords_plot(n2c,2)], ...
      [coords_plot(n1c,3), coords_plot(n2c,3)], ...
      'r-', 'LineWidth', 4);

% -----------------------------
% Clean appearance
% -----------------------------
axis off;
box off;
set(gcf, 'Color', 'w');

title(sprintf('DepRate = %.2f, Critical Bar = %d', DepRate, maxDCR_idx));
hold off;



% -----------------------------------------------------------------------
%  Helper function
% -----------------------------------------------------------------------
function s = sel_str(cond, s_true, s_false)
    if cond
        s = s_true;
    else
        s = s_false;
    end
end