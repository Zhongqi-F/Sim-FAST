clear all
close all
clc
tic

% =========================================================================
%  KIRIGAMI PEDESTRIAN BRIDGE - STRENGTH CHECK UNDER SELF-WEIGHT (SI Units)
%
%  Description : Structural analysis and LRFD strength check for a
%                Kirigami pedestrian bridge under dead load
%                only. This models the deployment scenario in which the
%                bridge must safely support its own self-weight.
%
%  Load:
%                DC : member self-weight (bar elements only)
%
%  Strength combination (AASHTO LRFD Table 3.4.1-1):
%                Strength I (DC only) : 1.25 * DC
%
%  Member strength checks (AASHTO LRFD Section 6):
%                Resistance factors   : Article 6.5.4.2
%                Tension resistance   : Article 6.8.2.1
%                Tension slenderness  : Article 6.8.4
%                Compression resist.  : Article 6.9.4.1
%                Compression slender. : Article 6.9.3
%
%  Purpose : Comparative structural performance study across bridge types.
%            Checks are limited to those necessary for fair cross-type
%            comparison under identical loading conditions.
%
%  Units   : SI throughout (N, m, Pa)
% =========================================================================


%% -------------------------------------------------------------------------
%  SECTION 1: GEOMETRY PARAMETERS
% -------------------------------------------------------------------------

L   = 2;    % Length of one repeating truss module (m)
gap = 0;    % Node gap offset within each module (m); 0 = no offset
N   = 8;    % Number of repeating modules along bridge span


%% -------------------------------------------------------------------------
%  SECTION 2: CROSS-SECTION PROPERTIES  (HSS 4x3x5/16, A500 Grade C)
%             Source: AISC Steel Construction Manual
% -------------------------------------------------------------------------

barE  = 200e9;      % Elastic modulus, steel (Pa)
Fy    = 345e6;      % Yield strength, A500 Gr C (Pa)
Fu    = 427e6;      % Tensile strength, A500 Gr C (Pa) — for 6.8.2.1-2

% AISC tabulated values (imperial)
t_in   = 0.291;     % Wall thickness (in)
A_in2  = 3.52;      % Cross-sectional area (in^2)
rx_in  = 1.42;      % Radius of gyration, strong axis x-x (in)
ry_in  = 1.13;      % Radius of gyration, weak axis y-y (in) <- governs
bt     = 7.31;      % Width-to-thickness ratio b/t
ht     = 10.7;      % Height-to-thickness ratio h/t

% Unit conversions: imperial -> SI
in2_to_m2 = (0.0254)^2;
in_to_m   = 0.0254;

barA = A_in2 * in2_to_m2;   % Cross-sectional area (m^2)
rx   = rx_in * in_to_m;     % Strong-axis radius of gyration (m)
ry   = ry_in * in_to_m;     % Weak-axis radius of gyration (m) <- used for KL/r

% Local slenderness check (AISC 360 Table B4.1a)
HSS = Check_HSS_Rect_Slenderness(bt, ht, barE, Fy, 'Forming', 'cold-formed');
fprintf("HSS local slenderness: b/t=%.2f, h/t=%.2f, lambda_r=%.2f -> %s\n", ...
        HSS.bt, HSS.ht, HSS.lambda_r, HSS.classStr);


%% -------------------------------------------------------------------------
%  SECTION 3: PANEL (SHEET) PROPERTIES
%             Panels are intentionally soft so that the truss carries all
%             global load; panel stiffness is set to ~1/1000 of steel.
% -------------------------------------------------------------------------

panel_E = 2e8;      % Panel elastic modulus (Pa)
panel_t = 0.01;     % Panel thickness (m)
panel_v = 0.3;      % Panel Poisson's ratio


%% -------------------------------------------------------------------------
%  SECTION 4: ASSEMBLY INITIALISATION
% -------------------------------------------------------------------------

assembly    = Assembly_Kirigami_Truss;
node        = Elements_Nodes;
cst         = Vec_Elements_CST;
rot_spr_4N  = Vec_Elements_RotSprings_4N;
rot_spr_3N  = CD_Elements_RotSprings_3N;
bar         = Vec_Elements_Bars;

assembly.cst        = cst;
assembly.node       = node;
assembly.bar        = bar;
assembly.rot_spr_4N = rot_spr_4N;
assembly.rot_spr_3N = rot_spr_3N;


%% -------------------------------------------------------------------------
%  SECTION 5: NODE COORDINATES
% -------------------------------------------------------------------------

% Start-face boundary nodes (x = 0)
node.coordinates_mat = [node.coordinates_mat;
    0, 0, 0;
    0, L, 0;
    0, 0, L;
    0, L, L;];

% Interior and end-plane nodes for each module
for i = 1:N
    node.coordinates_mat = [node.coordinates_mat;
        (L)*(i-1)+L/2, 0,   0;
        (L)*(i-1)+L/2, 0,   gap;
        (L)*(i-1)+L/2, L,   0;
        (L)*(i-1)+L/2, L,   gap;

        (L)*(i-1)+L/2, 0,   L;
        (L)*(i-1)+L/2, gap, L;
        (L)*(i-1)+L/2, L,   L;
        (L)*(i-1)+L/2, L,   L-gap;

        (L)*(i-1)+L/2, L/2, 0;
        (L)*(i-1)+L/2, L/2, L;
        (L)*(i-1)+L/2, 0,   L/2;
        (L)*(i-1)+L/2, L,   L/2;

        (L)*(i-1)+L,   0,   0;
        (L)*(i-1)+L,   L,   0;
        (L)*(i-1)+L,   0,   L;
        (L)*(i-1)+L,   L,   L;];
end


%% -------------------------------------------------------------------------
%  SECTION 6: PLOTTING SETUP
% -------------------------------------------------------------------------

plots = Plot_Kirigami_Truss;
plots.assembly     = assembly;
plots.displayRange = [-1; L*(N+1); -1; 3; -1; 3];
plots.viewAngle1   = 20;
plots.viewAngle2   = 20;

plots.Plot_Shape_Node_Number;


%% -------------------------------------------------------------------------
%  SECTION 7: CST PANEL ELEMENTS
% -------------------------------------------------------------------------

for i = 1:N
    cst.node_ijk_mat = [cst.node_ijk_mat;
        16*(i-1)+1,  16*(i-1)+5,  16*(i-1)+13;
        16*(i-1)+1,  16*(i-1)+13, 16*(i-1)+2;
        16*(i-1)+2,  16*(i-1)+7,  16*(i-1)+13;
        16*(i-1)+7,  16*(i-1)+13, 16*(i-1)+18;
        16*(i-1)+13, 16*(i-1)+17, 16*(i-1)+18;
        16*(i-1)+13, 16*(i-1)+5,  16*(i-1)+17;
        ];
end

cstNum    = size(cst.node_ijk_mat, 1);
cst.t_vec = panel_t * ones(cstNum, 1);
cst.E_vec = panel_E * ones(cstNum, 1);
cst.v_vec = panel_v * ones(cstNum, 1);

plots.Plot_Shape_CST_Number;


%% -------------------------------------------------------------------------
%  SECTION 8: BAR (TRUSS) ELEMENTS
% -------------------------------------------------------------------------

for i = 1:N
    bar.node_ij_mat = [bar.node_ij_mat;
        16*(i-1)+1,  16*(i-1)+2;
        16*(i-1)+1,  16*(i-1)+3;
        16*(i-1)+2,  16*(i-1)+4;
        16*(i-1)+3,  16*(i-1)+4;

        16*(i-1)+1,  16*(i-1)+5;
        16*(i-1)+5,  16*(i-1)+17;
        16*(i-1)+2,  16*(i-1)+7;
        16*(i-1)+7,  16*(i-1)+18;

        16*(i-1)+3,  16*(i-1)+9;
        16*(i-1)+9,  16*(i-1)+19;
        16*(i-1)+4,  16*(i-1)+11;
        16*(i-1)+11, 16*(i-1)+20;

        16*(i-1)+4,  16*(i-1)+12;
        16*(i-1)+12, 16*(i-1)+20;
        16*(i-1)+3,  16*(i-1)+10;
        16*(i-1)+10, 16*(i-1)+19;

        16*(i-1)+1,  16*(i-1)+6;
        16*(i-1)+6,  16*(i-1)+17;
        16*(i-1)+2,  16*(i-1)+8;
        16*(i-1)+8,  16*(i-1)+18;

        16*(i-1)+1,  16*(i-1)+15;
        16*(i-1)+3,  16*(i-1)+15;
        16*(i-1)+15, 16*(i-1)+17;
        16*(i-1)+15, 16*(i-1)+19;

        16*(i-1)+2,  16*(i-1)+16;
        16*(i-1)+16, 16*(i-1)+20;
        16*(i-1)+4,  16*(i-1)+16;
        16*(i-1)+16, 16*(i-1)+18;

        16*(i-1)+3,  16*(i-1)+14;
        16*(i-1)+4,  16*(i-1)+14;
        16*(i-1)+14, 16*(i-1)+20;
        16*(i-1)+14, 16*(i-1)+19;

        16*(i-1)+1,  16*(i-1)+13;
        16*(i-1)+2,  16*(i-1)+13;
        16*(i-1)+13, 16*(i-1)+18;
        16*(i-1)+13, 16*(i-1)+17;

        16*(i-1)+10, 16*(i-1)+14;
        16*(i-1)+11, 16*(i-1)+14;
        16*(i-1)+12, 16*(i-1)+16;
        16*(i-1)+8,  16*(i-1)+16;

        16*(i-1)+7,  16*(i-1)+13;
        16*(i-1)+5,  16*(i-1)+13;
        16*(i-1)+6,  16*(i-1)+15;
        16*(i-1)+15, 16*(i-1)+9;
        ];
end

% Closing bars on the terminal end face (x = N*L)
i = N + 1;
bar.node_ij_mat = [bar.node_ij_mat;
    16*(i-1)+1, 16*(i-1)+2;
    16*(i-1)+1, 16*(i-1)+3;
    16*(i-1)+2, 16*(i-1)+4;
    16*(i-1)+3, 16*(i-1)+4;];

barNum    = size(bar.node_ij_mat, 1);
bar.A_vec = barA * ones(barNum, 1);
bar.E_vec = barE * ones(barNum, 1);

plots.Plot_Shape_Bar_Number();
plots.Plot_Shape_Node_Number();


%% -------------------------------------------------------------------------
%  SECTION 9: 4-NODE ROTATIONAL SPRINGS
% -------------------------------------------------------------------------

for i = 1:N
    rot_spr_4N.node_ijkl_mat = [rot_spr_4N.node_ijkl_mat;
        16*(i-1)+1,  16*(i-1)+6,  16*(i-1)+15, 16*(i-1)+17;
        16*(i-1)+3,  16*(i-1)+9,  16*(i-1)+15, 16*(i-1)+19;

        16*(i-1)+1,  16*(i-1)+3,  16*(i-1)+15, 16*(i-1)+9;
        16*(i-1)+3,  16*(i-1)+1,  16*(i-1)+15, 16*(i-1)+6;
        16*(i-1)+6,  16*(i-1)+15, 16*(i-1)+17, 16*(i-1)+19;
        16*(i-1)+9,  16*(i-1)+15, 16*(i-1)+19, 16*(i-1)+17;

        16*(i-1)+3,  16*(i-1)+10, 16*(i-1)+14, 16*(i-1)+19;
        16*(i-1)+4,  16*(i-1)+11, 16*(i-1)+14, 16*(i-1)+20;

        16*(i-1)+11, 16*(i-1)+14, 16*(i-1)+4,  16*(i-1)+3;
        16*(i-1)+4,  16*(i-1)+14, 16*(i-1)+3,  16*(i-1)+10;
        16*(i-1)+10, 16*(i-1)+14, 16*(i-1)+19, 16*(i-1)+20;
        16*(i-1)+11, 16*(i-1)+14, 16*(i-1)+20, 16*(i-1)+19;

        16*(i-1)+2,  16*(i-1)+8,  16*(i-1)+16, 16*(i-1)+18;
        16*(i-1)+4,  16*(i-1)+12, 16*(i-1)+16, 16*(i-1)+20;

        16*(i-1)+2,  16*(i-1)+16, 16*(i-1)+4,  16*(i-1)+12;
        16*(i-1)+4,  16*(i-1)+16, 16*(i-1)+2,  16*(i-1)+8;
        16*(i-1)+8,  16*(i-1)+16, 16*(i-1)+18, 16*(i-1)+20;
        16*(i-1)+18, 16*(i-1)+16, 16*(i-1)+20, 16*(i-1)+12;

        16*(i-1)+2,  16*(i-1)+7,  16*(i-1)+13, 16*(i-1)+18;
        16*(i-1)+1,  16*(i-1)+5,  16*(i-1)+13, 16*(i-1)+17;

        16*(i-1)+1,  16*(i-1)+13, 16*(i-1)+2,  16*(i-1)+7;
        16*(i-1)+2,  16*(i-1)+13, 16*(i-1)+1,  16*(i-1)+5;
        16*(i-1)+5,  16*(i-1)+13, 16*(i-1)+17, 16*(i-1)+18;
        16*(i-1)+7,  16*(i-1)+13, 16*(i-1)+18, 16*(i-1)+17;
        ];
end

rotNum = size(rot_spr_4N.node_ijkl_mat, 1);
rot_spr_4N.rot_spr_K_vec = 1e8 * ones(rotNum, 1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_Spr_Number;


%% -------------------------------------------------------------------------
%  SECTION 10: 3-NODE ROTATIONAL SPRINGS
% -------------------------------------------------------------------------

for i = 1:N+1
    rot_spr_3N.node_ijk_mat = [rot_spr_3N.node_ijk_mat;
        16*(i-1)+1, 16*(i-1)+2, 16*(i-1)+4;
        16*(i-1)+2, 16*(i-1)+4, 16*(i-1)+3;
        16*(i-1)+4, 16*(i-1)+3, 16*(i-1)+1;
        16*(i-1)+3, 16*(i-1)+1, 16*(i-1)+2;
        ];
end

rot3Num = size(rot_spr_3N.node_ijk_mat, 1);
rot_spr_3N.rot_spr_K_vec = 1e8 * ones(rot3Num, 1);

plots.Plot_Shape_RotSpr_3N_Number();


%% -------------------------------------------------------------------------
%  SECTION 11: ASSEMBLY INITIALISATION
% -------------------------------------------------------------------------

assembly.Initialize_Assembly;


%% =========================================================================
%  LOAD DEFINITION
% =========================================================================

rho_steel = 7850;   % Steel density (kg/m^3)
g         = 9.81;   % Gravitational acceleration (m/s^2)

% DC: bar self-weight distributed to end nodes of each bar
LoadCase.DC = Build_LoadCase_DC_Bars(bar.node_ij_mat, node.coordinates_mat, ...
                                     barA, rho_steel, g);

Wbar = -sum(LoadCase.DC(:,4));
fprintf('DC total (bar self-weight) = %.2f N\n', Wbar);


%% =========================================================================
%  BOUNDARY CONDITIONS
%
%  Start face (x = 0) : fully pinned — Ux, Uy, Uz all fixed
%  End face   (x = NL): roller      — Ux free, Uy and Uz fixed
% =========================================================================

nr          = Solver_NR_Loading;
nr.assembly = assembly;
nodeNum     = size(node.coordinates_mat, 1);

% [nodeID, Ux, Uy, Uz]  (1 = fixed, 0 = free)
nr.supp = [
    1,       1, 1, 1;    % start face, bottom node 1
    2,       1, 1, 1;    % start face, bottom node 2
    16*N+1,  0, 1, 1;    % end face,   bottom node 1
    16*N+2,  0, 1, 1;    % end face,   bottom node 2
    ];


%% =========================================================================
%  MEMBER PARAMETERS
%
%  Effective length factor K = 1.0 (pin-pin, AASHTO 4.6.2.5).
%  Weak-axis ry governs buckling (ry < rx for HSS 4x3).
% =========================================================================

barNum_global = size(bar.node_ij_mat, 1);

Keff_vec      = ones(barNum_global, 1);
L0_vec_global = bar.L0_vec(:);
KL_vec_global = Keff_vec .* L0_vec_global;
r_vec_global  = ry * ones(barNum_global, 1);


%% =========================================================================
%  MEMBER TYPE CLASSIFICATION
%  TopChord    : both nodes at z = zMax
%  BottomChord : both nodes at z = zMin
%  Web         : all other members
% =========================================================================

Z    = node.coordinates_mat(:, 3);
zMax = max(Z);
zMin = min(Z);
zTol = max(1e-9, 1e-6 * max(1, abs(zMax - zMin)));

isTopNode    = abs(Z - zMax) <= zTol;
isBottomNode = abs(Z - zMin) <= zTol;

memberType_global = strings(barNum_global, 1);

for e = 1:barNum_global
    n1e = bar.node_ij_mat(e, 1);
    n2e = bar.node_ij_mat(e, 2);
    if     isTopNode(n1e)    && isTopNode(n2e)
        memberType_global(e) = "TopChord";
    elseif isBottomNode(n1e) && isBottomNode(n2e)
        memberType_global(e) = "BottomChord";
    else
        memberType_global(e) = "Web";
    end
end

fprintf("Member classification: TopChord=%d, BottomChord=%d, Web=%d\n", ...
    sum(memberType_global == "TopChord"),  ...
    sum(memberType_global == "BottomChord"), ...
    sum(memberType_global == "Web"));


%% =========================================================================
%  RESISTANCE FACTORS  (AASHTO LRFD 6.5.4.2)
%  Net section parameters assume welded connection (An = Ag, U = 1.0).
% =========================================================================

phi_ty  = 0.95;   % tension yielding in gross section    (6.5.4.2)
phi_cb  = 0.95;   % axial compression, steel only        (6.5.4.2)
phi_uf  = 0.80;   % tension fracture in net section      (6.5.4.2)

An_frac = barA;   % An = Ag: welded connection, no bolt holes
Rp_frac = 1.0;    % drilled/reamed holes                 (6.8.2.1)
U_shlag = 1.0;    % force transmitted to all elements    (6.8.2.1)


%% =========================================================================
%  STRENGTH CHECK  —  1.25 * DC
%
%  Single load combination representing the deployment scenario:
%  the bridge supports its own self-weight with the Strength I dead load
%  factor of 1.25 (AASHTO LRFD Table 3.4.1-1, DC maximum).
%
%  Member checks per AASHTO LRFD Section 6:
%    Tension     : 6.8.2.1  (gross yielding and net section fracture)
%    Compression : 6.9.4.1  (inelastic / elastic column curve)
%    Slenderness : 6.8.4 (tension), 6.9.3 (compression)
% =========================================================================

fprintf("\n=== Strength Check: 1.25 * DC ===\n");

nr.increStep = 1;
nr.iterMax   = 50;
nr.tol       = 1e-5;

nr.load = Build_Combo_Load(LoadCase, struct('DC', 1.25, 'PL', 0.0), nodeNum);
total_F = -sum(nr.load(:, 4));
fprintf("  Total factored vertical load (1.25 DC) = %.2f N\n", total_F);

Uhis  = nr.Solve;
U_end = squeeze(Uhis(end, :, :));
Uz    = U_end(:, 3);
maxDisp = -min(Uz);
fprintf("  Max downward displacement              = %.6f m\n", maxDisp);

truss_strain   = bar.Solve_Strain(node, U_end);
internal_force = truss_strain .* bar.E_vec .* bar.A_vec;

Pn      = NaN(barNum_global, 1);
phi     = NaN(barNum_global, 1);
phiPn   = NaN(barNum_global, 1);
DCR     = NaN(barNum_global, 1);
modeStr = cell(barNum_global, 1);
passYN  = false(barNum_global, 1);

for k = 1:barNum_global
    [passi, modeStri, Pni, phii, phiPni, DCRi] = Check_Truss_LRFD( ...
        internal_force(k), bar.A_vec(k), bar.E_vec(k), ...
        KL_vec_global(k),  r_vec_global(k), Fy, ...
        phi_ty, phi_cb, ...
        Fu, An_frac, Rp_frac, U_shlag, phi_uf, ...
        memberType_global(k));

    passYN(k)  = passi;
    modeStr{k} = modeStri;
    Pn(k)      = Pni;
    phi(k)     = phii;
    phiPn(k)   = phiPni;
    DCR(k)     = DCRi;
end

[maxDCR, kcrit] = max(DCR, [], 'omitnan');

fprintf("  Max D/C ratio = %.3f at bar #%d\n", maxDCR, kcrit);
fprintf("  Failure mode  = %s\n", modeStr{kcrit});

n1       = bar.node_ij_mat(kcrit, 1);
n2       = bar.node_ij_mat(kcrit, 2);
p1       = node.coordinates_mat(n1, :);
p2       = node.coordinates_mat(n2, :);
Lcrit    = norm(p2 - p1);
critType = memberType_global(kcrit);

fprintf("  Critical member #%d (%s): nodes %d-%d, L = %.3f m\n", ...
        kcrit, critType, n1, n2, Lcrit);
fprintf("    Pu = %.1f N  |  Pn = %.1f N  |  phi = %.2f  |  phi*Pn = %.1f N\n", ...
        internal_force(kcrit), Pn(kcrit), phi(kcrit), phiPn(kcrit));

% Slenderness report for critical member
KLr_crit   = KL_vec_global(kcrit) / r_vec_global(kcrit);
isPrimCrit = (critType == "TopChord" || critType == "BottomChord");
if internal_force(kcrit) < 0
    KLr_lim = 120 * isPrimCrit + 140 * (~isPrimCrit);
    fprintf("    KL/r = %.2f (AASHTO 6.9.3 limit = %d for %s)\n", ...
            KLr_crit, KLr_lim, critType);
    Po_crit    = Fy * bar.A_vec(kcrit);
    Pe_crit    = (pi^2 * barE * bar.A_vec(kcrit)) / KLr_crit^2;
    ratio_crit = Po_crit / Pe_crit;
    branchStr  = "inelastic [Eq. 6.9.4.1.1-1]";
    if ratio_crit > 2.25, branchStr = "elastic [Eq. 6.9.4.1.1-2]"; end
    fprintf("    Po/Pe = %.3f -> %s\n", ratio_crit, branchStr);
else
    Lr_lim = 200 * isPrimCrit + 240 * (~isPrimCrit);
    fprintf("    l/r  = %.2f (AASHTO 6.8.4 limit = %d for %s)\n", ...
            KLr_crit, Lr_lim, critType);
end

% Store results
Results.Strength_DC.U_end    = U_end;
Results.Strength_DC.total_F  = total_F;
Results.Strength_DC.maxDisp  = maxDisp;
Results.Strength_DC.maxDCR   = maxDCR;
Results.Strength_DC.critBar  = kcrit;
Results.Strength_DC.critType = critType;
Results.Strength_DC.critMode = modeStr{kcrit};
Results.Strength_DC.critNodes= [n1, n2];
Results.Strength_DC.critLen  = Lcrit;
Results.Strength_DC.passYN   = passYN;
Results.Strength_DC.DCR      = DCR;
Results.Strength_DC.Pu       = internal_force;
Results.Strength_DC.Pn       = Pn;
Results.Strength_DC.phi      = phi;
Results.Strength_DC.phiPn    = phiPn;
Results.Strength_DC.modeStr  = modeStr;


%% =========================================================================
%  CAPACITY SCAN  —  DC load multiplier beta
%
%  Purpose : Find the maximum DC load multiplier beta* at which the first
%            member reaches its LRFD design capacity (DCR = 1.0).
%            This defines the structural capacity under self-weight loading.
%
%  Method  : Apply beta * DC and increase beta until max(DCR) >= 1.0.
%            Two-stage scan: coarse (step 0.5) then fine (step 0.05).
%
%  Member checks identical to Strength Check above (AASHTO LRFD Section 6):
%            Tension     : 6.8.2.1
%            Compression : 6.9.4.1
%            phi values  : 6.5.4.2
%
%  Capacity/Weight = beta* × Wbar / Wbar = beta*
% =========================================================================

fprintf("\n=== Capacity Scan: beta * DC ===\n");

betaCoarse = 0.5:0.5:200;   % coarse scan range (adjust upper limit if needed)
coarseStep = 0.5;
fineStep   = 0.05;

Cap_DC           = struct();
Cap_DC.beta      = NaN;
Cap_DC.total_F   = NaN;
Cap_DC.maxDCR    = NaN;
Cap_DC.critBar   = NaN;
Cap_DC.critType  = "";
Cap_DC.critMode  = "";
Cap_DC.critNodes = [NaN, NaN];
Cap_DC.critLen   = NaN;
Cap_DC.converged = false;

hit        = false;
bracket_lo = NaN;
bracket_hi = NaN;

% -------------------------------------------------------------------------
% Stage 1: Coarse scan
% -------------------------------------------------------------------------
for a = 1:numel(betaCoarse)

    beta = betaCoarse(a);

    nr.load      = Build_Combo_Load(LoadCase, struct('DC', beta, 'PL', 0.0), nodeNum);
    total_F_scan = -sum(nr.load(:, 4));

    nr.increStep = 1;
    nr.iterMax   = 50;
    nr.tol       = 1e-5;

    try
        Uhis_scan = nr.Solve;
    catch ME
        fprintf("  beta=%.2f solver failed: %s\n", beta, ME.message);
        break
    end

    U_end_scan    = squeeze(Uhis_scan(end, :, :));
    truss_s       = bar.Solve_Strain(node, U_end_scan);
    iforce_scan   = truss_s .* bar.E_vec .* bar.A_vec;

    DCR_scan      = NaN(barNum_global, 1);
    modeStr_scan  = cell(barNum_global, 1);

    for k = 1:barNum_global
        [~, modeStr_scan{k}, ~, ~, ~, DCR_scan(k)] = Check_Truss_LRFD( ...
            iforce_scan(k), bar.A_vec(k), bar.E_vec(k), ...
            KL_vec_global(k), r_vec_global(k), Fy, ...
            phi_ty, phi_cb, ...
            Fu, An_frac, Rp_frac, U_shlag, phi_uf, ...
            memberType_global(k));
    end

    [maxDCR_scan, kcrit_scan] = max(DCR_scan, [], 'omitnan');

    fprintf("  beta=%6.2f  TotalF=%10.1f N  MaxD/C=%.3f  (bar#%d, %s)\n", ...
            beta, total_F_scan, maxDCR_scan, kcrit_scan, ...
            memberType_global(kcrit_scan));

    if maxDCR_scan >= 1.0
        hit        = true;
        bracket_hi = beta;
        bracket_lo = beta - coarseStep;
        fprintf("  >>> Bracket found (coarse): [%.2f, %.2f]\n", ...
                bracket_lo, bracket_hi);
        break
    end
end

% -------------------------------------------------------------------------
% Stage 2: Fine refinement
% -------------------------------------------------------------------------
if hit
    fprintf("  --- Refinement: [%.2f, %.2f], step=%.2f ---\n", ...
            bracket_lo, bracket_hi, fineStep);

    betaFine = bracket_lo:fineStep:bracket_hi;

    for a = 1:numel(betaFine)

        beta = betaFine(a);

        nr.load      = Build_Combo_Load(LoadCase, struct('DC', beta, 'PL', 0.0), nodeNum);
        total_F_scan = -sum(nr.load(:, 4));

        nr.increStep = 1;
        nr.iterMax   = 50;
        nr.tol       = 1e-5;

        try
            Uhis_scan = nr.Solve;
        catch
            continue
        end

        U_end_scan   = squeeze(Uhis_scan(end, :, :));
        truss_s      = bar.Solve_Strain(node, U_end_scan);
        iforce_scan  = truss_s .* bar.E_vec .* bar.A_vec;

        DCR_scan     = NaN(barNum_global, 1);
        modeStr_scan = cell(barNum_global, 1);

        for k = 1:barNum_global
            [~, modeStr_scan{k}, ~, ~, ~, DCR_scan(k)] = Check_Truss_LRFD( ...
                iforce_scan(k), bar.A_vec(k), bar.E_vec(k), ...
                KL_vec_global(k), r_vec_global(k), Fy, ...
                phi_ty, phi_cb, ...
                Fu, An_frac, Rp_frac, U_shlag, phi_uf, ...
                memberType_global(k));
        end

        [maxDCR_scan, kcrit_scan] = max(DCR_scan, [], 'omitnan');

        fprintf("    beta=%6.2f  TotalF=%10.1f N  MaxD/C=%.3f  (bar#%d, %s)\n", ...
                beta, total_F_scan, maxDCR_scan, kcrit_scan, ...
                memberType_global(kcrit_scan));

        if maxDCR_scan >= 1.0
            n1_scan = bar.node_ij_mat(kcrit_scan, 1);
            n2_scan = bar.node_ij_mat(kcrit_scan, 2);
            p1_scan = node.coordinates_mat(n1_scan, :);
            p2_scan = node.coordinates_mat(n2_scan, :);

            Cap_DC.beta      = beta;
            Cap_DC.total_F   = total_F_scan;
            Cap_DC.maxDCR    = maxDCR_scan;
            Cap_DC.critBar   = kcrit_scan;
            Cap_DC.critType  = memberType_global(kcrit_scan);
            Cap_DC.critMode  = string(modeStr_scan{kcrit_scan});
            Cap_DC.critNodes = [n1_scan, n2_scan];
            Cap_DC.critLen   = norm(p2_scan - p1_scan);
            Cap_DC.converged = true;

            fprintf("  >>> Refined capacity: beta*=%.2f, TotalF=%.1f N, MaxD/C=%.3f\n", ...
                    beta, total_F_scan, maxDCR_scan);
            break
        end
    end
end

% -------------------------------------------------------------------------
% Report Capacity/Weight
% -------------------------------------------------------------------------
fprintf("\n=== CAPACITY / WEIGHT RESULT ===\n");

if Cap_DC.converged
    Cap_to_Weight = Cap_DC.beta;   % = beta* × Wbar / Wbar = beta*

    fprintf("  beta*             = %.2f\n",       Cap_DC.beta);
    fprintf("  Capacity (beta*×W)= %.1f N\n",     Cap_DC.total_F);
    fprintf("  Self-weight (W)   = %.1f N\n",     Wbar);
    fprintf("  Capacity/Weight   = %.3f\n",        Cap_to_Weight);
    fprintf("  Governing bar     = #%d (%s)\n",   Cap_DC.critBar, Cap_DC.critType);
    fprintf("  Failure mode      = %s\n",          Cap_DC.critMode);
else
    Cap_to_Weight = NaN;
    fprintf("  [WARN] Capacity not reached within scan range.\n");
    fprintf("         Increase upper limit of betaCoarse.\n");
end

fprintf("================================\n");


%% =========================================================================
%  PLOTS
% =========================================================================

sigma_plot = Results.Strength_DC.Pu ./ bar.A_vec(:);   % axial stress (Pa)

plots.Plot_Shape_Bar_Stress(sigma_plot);
plots.Plot_Shape_Bar_Failure(Results.Strength_DC.passYN);

Plot_Highlight_Controlled_Bar( ...
    node.coordinates_mat, bar.node_ij_mat, kcrit, ...
    'ShowAllBars', true, ...
    'HiColor',     [1 0 0], ...
    'HiLineWidth', 3.5, ...
    'View',        [plots.viewAngle1, plots.viewAngle2], ...
    'Title',       sprintf("Controlling member: Strength DC (bar #%d)", kcrit));


%% =========================================================================
%  SYSTEM PERFORMANCE SUMMARY
% =========================================================================

L_unit     = L;
span_total = L_unit * N;

K_system = total_F / maxDisp;   % system stiffness under 1.25 DC (N/m)
K_to_W   = K_system / Wbar;     % stiffness-to-weight ratio

fprintf("\n===== SYSTEM PERFORMANCE SUMMARY =====\n");
fprintf("  Bar self-weight (DC)           = %10.1f N\n",       Wbar);
fprintf("  Factored load (1.25 DC)        = %10.1f N\n",       total_F);
fprintf("  Max displacement               = %10.6f m\n",       maxDisp);
fprintf("  System stiffness K             = %10.2e N/m\n",     K_system);
fprintf("  Stiffness-to-weight ratio      = %10.4f (N/m)/N\n", K_to_W);
fprintf("  Max D/C ratio (Strength DC)    = %10.3f\n",         maxDCR);
fprintf("  Governing member               = bar #%d (%s)\n",   kcrit, critType);
fprintf("  All members pass (DCR<=1.0)    = %s\n",             mat2str(all(passYN)));
fprintf("  Capacity/Weight ratio          = %10.3f\n", Cap_to_Weight);
fprintf("=======================================\n");

toc