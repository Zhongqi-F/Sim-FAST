clear all
close all
clc
tic

% =========================================================================
%  KIRIGAMI PEDESTRIAN BRIDGE - AASHTO LRFD STRUCTURAL ANALYSIS (SI Units)
%
%  Description : Structural analysis + LRFD-style member utilization check
%                for a Kirigami-inspired pedestrian truss bridge under:
%                - DC: member self-weight
%                - PL: pedestrian live load (uniform area load)
%
%  Load combinations (AASHTO LRFD Table 3.4.1-1):
%                Service I   : 1.00 DC + 1.00 PL
%                Strength I  : 1.25 DC + 1.75 PL   (DC max)
%                              0.90 DC + 1.75 PL   (DC min)
%
%  Design/check references:
%                - Loads & factors      : AASHTO LRFD Section 3
%                - Resistance factors   : AASHTO LRFD 6.5.4.2
%                - Tension resistance   : AASHTO LRFD 6.8.2.1
%                - Tension slenderness  : AASHTO LRFD 6.8.4
%                - Compression resist.  : AASHTO LRFD 6.9.4.1
%                - Compression slender. : AASHTO LRFD 6.9.3
%
%  TODO (requires Section 2 text):
%                - Deflection limit: verify L/360 per AASHTO LRFD 2.5.2.6.2
%                - Dead load DC: deck panel and railing weights per AASHTO 3.5.1
%                  currently only bar self-weight is included
%
%  TODO (requires Guide Spec for Pedestrian Bridges Section 6):
%                - Vibration/dynamic serviceability check not implemented
%                  AASHTO Guide Spec requires vertical natural frequency >= 3.0 Hz
%                  or detailed dynamic analysis; needs eigenvalue analysis eig(K,M)
%
%  Units       : SI throughout (N, m, Pa)
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
Fy    = 345e6;      % Yield strength, A500 Gr C ~ Grade 50 (Pa)
Fu    = 427e6;      % Tensile strength, A500 Gr C (Pa)  — used in 6.8.2.1-2

% AISC tabulated values (imperial)
t_in   = 0.291;     % Wall thickness (in)
A_in2  = 3.52;      % Cross-sectional area (in^2)
rx_in  = 1.42;      % Radius of gyration, strong axis x-x (in)
ry_in  = 1.13;      % Radius of gyration, weak axis y-y (in)  <- governs buckling
bt     = 7.31;      % Width-to-thickness ratio b/t
ht     = 10.7;      % Height-to-thickness ratio h/t

% Unit conversions: imperial -> SI
in2_to_m2 = (0.0254)^2;
in_to_m   = 0.0254;

barA = A_in2 * in2_to_m2;   % Cross-sectional area (m^2)
rx   = rx_in * in_to_m;     % Strong-axis radius of gyration (m)
ry   = ry_in * in_to_m;     % Weak-axis radius of gyration (m)  <- used for KL/r

% Local slenderness check (AISC 360 Table B4.1a)
HSS = Check_HSS_Rect_Slenderness(bt, ht, barE, Fy, 'Forming', 'cold-formed');
fprintf("HSS local slenderness: b/t=%.2f, h/t=%.2f, lambda_r=%.2f -> %s\n", ...
        HSS.bt, HSS.ht, HSS.lambda_r, HSS.classStr);


%% -------------------------------------------------------------------------
%  SECTION 3: PANEL (SHEET) PROPERTIES
% -------------------------------------------------------------------------

panel_E = 2e8;      % Panel elastic modulus (Pa)  ~200 MPa (soft)
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
plots.assembly    = assembly;
plots.displayRange = [-1; L*(N+1); -1; 3; -1; 3];
plots.viewAngle1  = 20;
plots.viewAngle2  = 20;

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

cstNum        = size(cst.node_ijk_mat, 1);
cst.t_vec     = panel_t   * ones(cstNum, 1);
cst.E_vec     = panel_E   * ones(cstNum, 1);
cst.v_vec     = panel_v   * ones(cstNum, 1);

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

barNum      = size(bar.node_ij_mat, 1);
bar.A_vec   = barA  * ones(barNum, 1);
bar.E_vec   = barE  * ones(barNum, 1);

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


%% -------------------------------------------------------------------------
%  SECTION 12: MATERIAL AND ANALYSIS CONTROL FLAGS
% -------------------------------------------------------------------------

rho_steel = 7850;   % Steel density (kg/m^3)
g         = 9.81;   % Gravitational acceleration (m/s^2)

DO_CODE_CHECK = true;
DO_CAPACITY   = true;


%% =========================================================================
%  LOAD DEFINITION
%  AASHTO LRFD Section 3
% =========================================================================

%% -------------------------------------------------------------------------
%  STEP 1A: Identify deck (bottom) nodes for live load application
% -------------------------------------------------------------------------

Z_coords  = node.coordinates_mat(:, 3);
zMin_deck = min(Z_coords);
zTol_deck = max(1e-9, 1e-6 * max(1, abs(zMin_deck)));

bottomNodeIDs = find(abs(Z_coords - zMin_deck) <= zTol_deck);
bottomNodeIDs = unique(bottomNodeIDs);

fprintf("PL bottom nodes detected (z=zMin): %d nodes\n", numel(bottomNodeIDs));
fprintf("Bottom node IDs (z=zMin):\n");
disp(bottomNodeIDs(:).');


%% -------------------------------------------------------------------------
%  STEP 1B: Build individual load case vectors
% -------------------------------------------------------------------------

LoadCase = struct();

% DC: bar self-weight distributed to end nodes of each bar
% TODO (AASHTO 3.5.1 — requires Section 3 text):
%   Deck panel self-weight and railing/barrier weights should also be
%   included in DC once section properties are confirmed.
LoadCase.DC = Build_LoadCase_DC_Bars(bar.node_ij_mat, node.coordinates_mat, ...
                                     barA, rho_steel, g);

% PL: uniform pedestrian pressure (AASHTO 3.6.1.6: 3.6 kPa)
% TODO (geometry verification required):
%   width_deck = L assumed equal to module length.
%   Verify this equals effective load-bearing deck width for Kirigami geometry.
L_unit     = L;
span_total = L_unit * N;
width_deck = L_unit;

qPL = 3.6e3;   % Pa — AASHTO LRFD 3.6.1.6 pedestrian live load
LoadCase.PL = Build_LoadCase_PL_Uniform(bottomNodeIDs, qPL, span_total, width_deck);
fprintf("Check PL total (q*span*width) = %.2f N\n", qPL*span_total*width_deck);

Wbar_check = -sum(LoadCase.DC(:,4));
PL_check   = -sum(LoadCase.PL(:,4));
fprintf('DC total (bar self-weight)  = %.2f N\n', Wbar_check);
fprintf('PL total (alpha = 1.0)      = %.2f N\n', PL_check);


%% =========================================================================
%  STEP 2: LOAD COMBINATIONS  (AASHTO LRFD Table 3.4.1-1)
% =========================================================================

Combos = struct();

Combos.Service_1.name    = "Service_1";
Combos.Service_1.factors = struct('DC', 1.00, 'PL', 1.00);

Combos.Strength_1a.name    = "Strength_1a";
Combos.Strength_1a.factors = struct('DC', 1.25, 'PL', 1.75);   % DC maximum

Combos.Strength_1b.name    = "Strength_1b";
Combos.Strength_1b.factors = struct('DC', 0.90, 'PL', 1.75);   % DC minimum


%% =========================================================================
%  STEP 3: BOUNDARY CONDITIONS
%
%  Simple support at both abutments:
%    Start face (x=0): fully pinned — Ux, Uy, Uz all fixed
%    End face (x=N*L): sliding pin  — Uy, Uz fixed, Ux released (roller)
%
%  All four corner nodes at each end face are constrained to prevent
%  rigid-body motion and end-face distortion.
% =========================================================================

nr           = Solver_NR_Loading;
nr.assembly  = assembly;
nodeNum      = size(node.coordinates_mat, 1);

% [nodeID, Ux, Uy, Uz]  (1 = fixed, 0 = free)
nr.supp = [
    % Start face — fully pinned
    1,        1, 1, 1;
    2,        1, 1, 1;
    % End face — roller (Ux free, Uy and Uz fixed)
    16*N+1,   0, 1, 1;
    16*N+2,   0, 1, 1;
    ];


%% =========================================================================
%  STEP 4: MEMBER BUCKLING PARAMETERS
%
%  Effective length factor K = 1.0 (pin-pin assumption, AASHTO 4.6.2.5).
%  Weak-axis radius of gyration ry governs (ry < rx for HSS 4x3).
% =========================================================================

barNum_global = size(bar.node_ij_mat, 1);

Keff_vec      = ones(barNum_global, 1);          % K = 1.0 for all members
L0_vec_global = bar.L0_vec(:);                   % undeformed member lengths (m)
KL_vec_global = Keff_vec .* L0_vec_global;       % effective lengths KL (m)
r_vec_global  = ry * ones(barNum_global, 1);     % weak-axis ry governs (m)


%% =========================================================================
%  STEP 5: MEMBER TYPE CLASSIFICATION
%          TopChord    : both nodes at z = zMax
%          BottomChord : both nodes at z = zMin
%          Web         : all other members
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
%  STEP 6: RESISTANCE FACTORS AND NET SECTION PARAMETERS
%          AASHTO LRFD 6.5.4.2
% =========================================================================

% Resistance factors — AASHTO LRFD 6.5.4.2
phi_ty  = 0.95;   % tension yielding in gross section       (6.5.4.2)
phi_cb  = 0.95;   % axial compression, steel only           (6.5.4.2)
phi_uf  = 0.80;   % tension fracture in net section         (6.5.4.2)

% Net section fracture parameters — AASHTO LRFD 6.8.2.1-2
% Assumption: welded HSS connection, force transmitted to all elements
An_frac  = barA;   % An = Ag: welded connection, no bolt holes assumed
Rp_frac  = 1.0;    % drilled/reamed holes per 6.8.2.1
U_shlag  = 1.0;    % U = 1.0: force transmitted to all elements per 6.8.2.1


%% =========================================================================
%  PATH 1: PER-COMBINATION CODE CHECK
%          For each load combination:
%            1. Assemble factored nodal load vector
%            2. Solve for displacements (Newton-Raphson)
%            3. Recover bar forces
%            4. Compute LRFD demand-to-capacity ratio (DCR)
%               per AASHTO LRFD 6.8.2.1 (tension) and 6.9.4.1 (compression)
%            5. Report maximum DCR and governing member
%            6. Check KL/r slenderness limits per 6.9.3 and 6.8.4
% =========================================================================

Results = struct();

if DO_CODE_CHECK

    comboNames = fieldnames(Combos);

    for c = 1:numel(comboNames)

        comboField = comboNames{c};
        combo      = Combos.(comboField);

        fprintf("\n=== Running combination: %s ===\n", combo.name);

        nr.increStep = 1;
        nr.iterMax   = 50;
        nr.tol       = 1e-5;

        nr.load  = Build_Combo_Load(LoadCase, combo.factors, nodeNum);
        total_F  = -sum(nr.load(:, 4));
        fprintf("  Total factored vertical load = %.2f N\n", total_F);

        Uhis    = nr.Solve;
        U_end   = squeeze(Uhis(end, :, :));
        Uz      = U_end(:, 3);
        maxDisp = -min(Uz);

        Results.(comboField).U_end = U_end;

        fprintf("  Max downward displacement    = %.6f m\n", maxDisp);

        truss_strain   = bar.Solve_Strain(node, U_end);
        internal_force = truss_strain .* bar.E_vec .* bar.A_vec;

        barNum  = numel(internal_force);
        Pn      = NaN(barNum, 1);
        phi     = NaN(barNum, 1);
        phiPn   = NaN(barNum, 1);
        DCR     = NaN(barNum, 1);
        modeStr = cell(barNum, 1);
        passYN  = false(barNum, 1);

        for k = 1:barNum
            % AASHTO LRFD member check:
            %   Tension   : 6.8.2.1  (gross yielding and net fracture)
            %   Compression: 6.9.4.1 (inelastic/elastic column curve)
            %   phi values : 6.5.4.2
            %   Slenderness: 6.8.4 (tension), 6.9.3 (compression)
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

        if isnan(maxDCR) || isempty(kcrit)
            fprintf("  [WARN] Max D/C = NaN. Skipping critical member report.\n");
            Results.(comboField).total_F = total_F;
            Results.(comboField).maxDisp = maxDisp;
            Results.(comboField).maxDCR  = maxDCR;
            Results.(comboField).passYN  = passYN;
            continue
        end

        fprintf("  Max D/C ratio = %.3f at bar #%d, failure mode = %s\n", ...
                maxDCR, kcrit, modeStr{kcrit});

        n1    = bar.node_ij_mat(kcrit, 1);
        n2    = bar.node_ij_mat(kcrit, 2);
        p1    = node.coordinates_mat(n1, :);
        p2    = node.coordinates_mat(n2, :);
        Lcrit = norm(p2 - p1);
        critType = memberType_global(kcrit);

        fprintf("  Critical member #%d (%s): nodes %d-%d, L = %.3f m\n", ...
                kcrit, critType, n1, n2, Lcrit);
        fprintf("    Pu = %.1f N  |  Pn = %.1f N  |  phi = %.2f  |  phi*Pn = %.1f N\n", ...
                internal_force(kcrit), Pn(kcrit), phi(kcrit), phiPn(kcrit));

        % KL/r slenderness report for critical member
        KLr_crit  = KL_vec_global(kcrit) / r_vec_global(kcrit);
        isPrimCrit = (critType == "TopChord" || critType == "BottomChord");
        if internal_force(kcrit) < 0
            KLr_lim = 120 * isPrimCrit + 140 * (~isPrimCrit);
            fprintf("    KL/r = %.2f (AASHTO 6.9.3 limit = %d for %s)\n", ...
                    KLr_crit, KLr_lim, critType);
        else
            Lr_lim = 200 * isPrimCrit + 240 * (~isPrimCrit);
            fprintf("    l/r  = %.2f (AASHTO 6.8.4 limit = %d for %s)\n", ...
                    KLr_crit, Lr_lim, critType);
        end

        % Po/Pe slenderness indicator (AASHTO 6.9.4.1 nomenclature)
        if internal_force(kcrit) < 0
            Po_crit    = Fy * bar.A_vec(kcrit);
            Pe_crit    = (pi^2 * barE * bar.A_vec(kcrit)) / KLr_crit^2;
            ratio_crit = Po_crit / Pe_crit;
            if ratio_crit <= 2.25
                branchStr = "inelastic [Eq. 6.9.4.1.1-1]";
            else
                branchStr = "elastic   [Eq. 6.9.4.1.1-2]";
            end
            fprintf("    Po/Pe = %.3f -> %s\n", ratio_crit, branchStr);
        end

        Results.(comboField).total_F   = total_F;
        Results.(comboField).maxDisp   = maxDisp;
        Results.(comboField).maxDCR    = maxDCR;
        Results.(comboField).critBar   = kcrit;
        Results.(comboField).critMode  = modeStr{kcrit};
        Results.(comboField).critType  = critType;
        Results.(comboField).critNodes = [n1, n2];
        Results.(comboField).critLen   = Lcrit;
        Results.(comboField).passYN    = passYN;
        Results.(comboField).DCR       = DCR;
        Results.(comboField).Pu        = internal_force;
        Results.(comboField).Pn        = Pn;
        Results.(comboField).phi       = phi;
        Results.(comboField).phiPn     = phiPn;
        Results.(comboField).modeStr   = modeStr;

    end

    % Global controlling member across all combinations
    maxDCR_all = -Inf;
    govCombo   = "";
    govBar     = NaN;

    for c = 1:numel(comboNames)
        thisField = comboNames{c};
        if isfield(Results, thisField) && isfield(Results.(thisField), 'maxDCR')
            thisMax = Results.(thisField).maxDCR;
            if ~isnan(thisMax) && thisMax > maxDCR_all
                maxDCR_all = thisMax;
                govCombo   = string(thisField);
                govBar     = Results.(thisField).critBar;
            end
        end
    end

    fprintf("\n=== GLOBAL CONTROLLING MEMBER (PATH 1) ===\n");
    fprintf("  Governing combo = %s\n", govCombo);
    fprintf("  Governing bar   = #%d\n", govBar);
    fprintf("  Max D/C         = %.3f\n", maxDCR_all);
    fprintf("==========================================\n");

end


%% =========================================================================
%  PATH 2: LRFD CAPACITY SCAN  (Strength I, DC maximum and DC minimum)
%
%  Factored live load scaled by multiplier alpha:
%    DC_max case:  1.25*DC + alpha*(1.75*PL)
%    DC_min case:  0.90*DC + alpha*(1.75*PL)
%
%  Same phi values as PATH 1 (AASHTO 6.5.4.2) for consistency.
% =========================================================================

if DO_CAPACITY

    fprintf("\n=============================\n");
    fprintf("CAPACITY SCAN (LRFD): Strength I, DC_max and DC_min\n");
    fprintf("=============================\n");

    DC_cases  = [1.25, 0.90];
    DC_labels = ["DC_max(1.25)", "DC_min(0.90)"];

    alphaCoarse = 0.5:0.5:50;
    coarseStep  = 0.5;
    fineStep    = 0.05;

    Cap_all = struct();

    for dc_idx = 1:numel(DC_cases)

        DC_factor = DC_cases(dc_idx);
        DC_label  = DC_labels(dc_idx);

        fprintf("\n--- Scanning %s ---\n", DC_label);

        Cap            = struct();
        Cap.DC_factor  = DC_factor;
        Cap.alpha      = NaN;
        Cap.total_F    = NaN;
        Cap.maxDCR     = NaN;
        Cap.critBar    = NaN;
        Cap.critMode   = "";
        Cap.critType   = "";
        Cap.critNodes  = [NaN, NaN];
        Cap.critLen    = NaN;
        Cap.converged  = false;

        hit        = false;
        bracket_lo = NaN;
        bracket_hi = NaN;

        % --- Coarse scan ---
        for a = 1:numel(alphaCoarse)

            alpha = alphaCoarse(a);

            nr.load  = Build_Combo_Load(LoadCase, ...
                           struct('DC', DC_factor, 'PL', 1.75*alpha), nodeNum);
            total_F  = -sum(nr.load(:, 4));

            nr.increStep = 1;
            nr.iterMax   = 50;
            nr.tol       = 1e-5;

            try
                Uhis = nr.Solve;
            catch ME
                fprintf("  alpha=%.2f solver failed: %s\n", alpha, ME.message);
                break
            end

            U_end          = squeeze(Uhis(end, :, :));
            truss_strain   = bar.Solve_Strain(node, U_end);
            internal_force = truss_strain .* bar.E_vec .* bar.A_vec;

            barNum_  = numel(internal_force);
            DCR_     = NaN(barNum_, 1);
            modeStr_ = cell(barNum_, 1);

            for k = 1:barNum_
                [~, modeStr_{k}, ~, ~, ~, DCR_(k)] = Check_Truss_LRFD( ...
                    internal_force(k), bar.A_vec(k), bar.E_vec(k), ...
                    KL_vec_global(k),  r_vec_global(k), Fy, ...
                    phi_ty, phi_cb, ...
                    Fu, An_frac, Rp_frac, U_shlag, phi_uf, ...
                    memberType_global(k));
            end

            [maxDCR_, kcrit_] = max(DCR_, [], 'omitnan');
            critType_         = memberType_global(kcrit_);

            fprintf("  alpha=%5.2f  TotalF=%10.1f N  MaxD/C=%.3f  (bar#%d, %s, %s)\n", ...
                    alpha, total_F, maxDCR_, kcrit_, modeStr_{kcrit_}, critType_);

            if maxDCR_ >= 1.0
                hit        = true;
                bracket_hi = alpha;
                bracket_lo = alpha - coarseStep;
                fprintf("  >>> Capacity bracket found (coarse): [%.2f, %.2f] [%s]\n", ...
                        bracket_lo, bracket_hi, DC_label);
                break
            end
        end

        % --- Fine refinement ---
        if hit
            fprintf("  --- Refinement: [%.2f, %.2f], step = %.2f ---\n", ...
                    bracket_lo, bracket_hi, fineStep);

            alphaFine = bracket_lo:fineStep:bracket_hi;

            for a = 1:numel(alphaFine)

                alpha = alphaFine(a);

                nr.load  = Build_Combo_Load(LoadCase, ...
                               struct('DC', DC_factor, 'PL', 1.75*alpha), nodeNum);
                total_F  = -sum(nr.load(:, 4));

                nr.increStep = 1;
                nr.iterMax   = 50;
                nr.tol       = 1e-5;

                try
                    Uhis = nr.Solve;
                catch
                    continue
                end

                U_end          = squeeze(Uhis(end, :, :));
                truss_strain   = bar.Solve_Strain(node, U_end);
                internal_force = truss_strain .* bar.E_vec .* bar.A_vec;

                barNum_  = numel(internal_force);
                DCR_     = NaN(barNum_, 1);
                modeStr_ = cell(barNum_, 1);

                for k = 1:barNum_
                    [~, modeStr_{k}, ~, ~, ~, DCR_(k)] = Check_Truss_LRFD( ...
                        internal_force(k), bar.A_vec(k), bar.E_vec(k), ...
                        KL_vec_global(k),  r_vec_global(k), Fy, ...
                        phi_ty, phi_cb, ...
                        Fu, An_frac, Rp_frac, U_shlag, phi_uf, ...
                        memberType_global(k));
                end

                [maxDCR_, kcrit_] = max(DCR_, [], 'omitnan');

                fprintf("    alpha=%5.2f  TotalF=%10.1f N  MaxD/C=%.3f  (bar#%d, %s, %s)\n", ...
                        alpha, total_F, maxDCR_, kcrit_, modeStr_{kcrit_}, ...
                        memberType_global(kcrit_));

                if maxDCR_ >= 1.0
                    n1_   = bar.node_ij_mat(kcrit_, 1);
                    n2_   = bar.node_ij_mat(kcrit_, 2);
                    p1_   = node.coordinates_mat(n1_, :);
                    p2_   = node.coordinates_mat(n2_, :);

                    Cap.alpha     = alpha;
                    Cap.total_F   = total_F;
                    Cap.maxDCR    = maxDCR_;
                    Cap.critBar   = kcrit_;
                    Cap.critMode  = string(modeStr_{kcrit_});
                    Cap.critType  = memberType_global(kcrit_);
                    Cap.critNodes = [n1_, n2_];
                    Cap.critLen   = norm(p2_ - p1_);
                    Cap.converged = true;

                    fprintf("  >>> Refined capacity: alpha=%.2f, TotalF=%.1f N, MaxD/C=%.3f [%s]\n", ...
                            alpha, total_F, maxDCR_, DC_label);
                    break
                end
            end
        end

        if DC_factor == 1.25
            Cap_all.DC_max = Cap;
        else
            Cap_all.DC_min = Cap;
        end

    end

    % Final capacity comparison
    fprintf("\n=== FINAL CAPACITY RESULT (LRFD, Strength I) ===\n");

    for dc_idx = 1:numel(DC_cases)
        if DC_cases(dc_idx) == 1.25
            Cap_i = Cap_all.DC_max;
        else
            Cap_i = Cap_all.DC_min;
        end

        if Cap_i.converged
            fprintf("  [%s]  alpha* = %.2f  |  TotalF* = %.1f N  |  MaxD/C = %.3f\n", ...
                    DC_labels(dc_idx), Cap_i.alpha, Cap_i.total_F, Cap_i.maxDCR);
            fprintf("         Governing bar #%d (%s), mode = %s\n", ...
                    Cap_i.critBar, Cap_i.critType, Cap_i.critMode);
        else
            fprintf("  [%s]  Did not converge or capacity not reached in scan range.\n", ...
                    DC_labels(dc_idx));
        end
    end

    if Cap_all.DC_max.converged && Cap_all.DC_min.converged
        if Cap_all.DC_max.total_F <= Cap_all.DC_min.total_F
            Cap = Cap_all.DC_max;
            fprintf("\n>>> GOVERNING: DC_max (1.25) controls  |  TotalF* = %.1f N\n", Cap.total_F);
        else
            Cap = Cap_all.DC_min;
            fprintf("\n>>> GOVERNING: DC_min (0.90) controls  |  TotalF* = %.1f N\n", Cap.total_F);
        end
    elseif Cap_all.DC_max.converged
        Cap = Cap_all.DC_max;
        fprintf("\n>>> Only DC_max converged. Using DC_max result.\n");
    elseif Cap_all.DC_min.converged
        Cap = Cap_all.DC_min;
        fprintf("\n>>> Only DC_min converged. Using DC_min result.\n");
    else
        Cap.converged = false;
        fprintf("\n>>> [WARN] Neither DC case reached capacity within scan range.\n");
    end

end


%% =========================================================================
%  POST-PROCESSING: PLOTS
% =========================================================================

PLOT_CASE = "Strength_1a";

if isfield(Results, PLOT_CASE)

    DCR_plot   = Results.(PLOT_CASE).DCR;
    pass_plot  = Results.(PLOT_CASE).passYN;
    Pu_plot    = Results.(PLOT_CASE).Pu;
    A_plot     = bar.A_vec(:);
    sigma_plot = Pu_plot ./ A_plot;

    [maxDCR_plot, kcrit_plot] = max(DCR_plot, [], 'omitnan');

    fprintf("\n=== PLOTTING CASE: %s ===\n", PLOT_CASE);
    fprintf("  Max D/C = %.3f at bar #%d (%s)\n", ...
            maxDCR_plot, kcrit_plot, Results.(PLOT_CASE).critType);

    % Po/Pe slenderness indicator for plot case critical member (AASHTO 6.9.4.1)
    KLr_crit  = KL_vec_global(kcrit_plot) / r_vec_global(kcrit_plot);
    Po_plot   = Fy * bar.A_vec(kcrit_plot);
    Pe_plot   = (pi^2 * barE * bar.A_vec(kcrit_plot)) / KLr_crit^2;
    ratio_plot = Po_plot / Pe_plot;

    if ratio_plot <= 2.25
        branchStr = "inelastic [Eq. 6.9.4.1.1-1]";
    else
        branchStr = "elastic   [Eq. 6.9.4.1.1-2]";
    end
    fprintf("  Governing slenderness: KL/r = %.2f,  Po/Pe = %.3f  (%s)\n", ...
            KLr_crit, ratio_plot, branchStr);

    plots.Plot_Shape_Bar_Stress(sigma_plot);
    plots.Plot_Shape_Bar_Failure(pass_plot);

    if isfield(Results, PLOT_CASE) && isfield(Results.(PLOT_CASE), 'critBar') ...
            && ~isempty(Results.(PLOT_CASE).critBar) && ~isnan(Results.(PLOT_CASE).critBar)

        barID = Results.(PLOT_CASE).critBar;

        Plot_Highlight_Controlled_Bar( ...
            node.coordinates_mat, bar.node_ij_mat, barID, ...
            'ShowAllBars',  true, ...
            'HiColor',      [1 0 0], ...
            'HiLineWidth',  3.5, ...
            'View',         [plots.viewAngle1, plots.viewAngle2], ...
            'Title',        sprintf("Controlling member: %s (bar #%d)", PLOT_CASE, barID));
    else
        fprintf("[WARN] No valid critBar for PLOT_CASE=%s. Skip highlight.\n", PLOT_CASE);
    end

else
    fprintf("\n[WARN] Results does not contain '%s'. Skipping plots.\n", PLOT_CASE);
end


%% =========================================================================
%  SYSTEM-LEVEL PERFORMANCE METRICS
% =========================================================================

L_unit     = L;
span_total = L_unit * N;
width_deck = L_unit;   % TODO: verify effective deck width for Kirigami geometry

if abs(span_total - L_unit) < 1e-12
    warning("span_total equals L_unit. Check N and naming.");
end

Wbar = -sum(LoadCase.DC(:,4));

delta_service = Results.Service_1.maxDisp;
F_service     = Results.Service_1.total_F;
K_service     = F_service / delta_service;
K_to_W        = K_service / Wbar;

% =========================================================================
% DEFLECTION CHECK
% TODO (AASHTO LRFD 2.5.2.6.2 — requires Section 2 text):
%   Deflection limit L/360 used as placeholder for pedestrian live load.
%   Confirm whether L/360 or another limit applies per AASHTO 2.5.2.6.2
%   and/or the AASHTO Guide Spec for Pedestrian Bridges.
%   The deflection check uses TOTAL SPAN (span_total), not module length.
% =========================================================================

nr_LL           = nr;
nr_LL.increStep = 1;
nr_LL.iterMax   = 50;
nr_LL.tol       = 1e-5;

% LL-only: DC=0, PL=1.0 (unfactored)
% TODO: Confirm this load case per AASHTO 2.5.2.6.2
nr_LL.load = Build_Combo_Load(LoadCase, struct('DC', 0.0, 'PL', 1.0), nodeNum);
Uhis_LL    = nr_LL.Solve;
U_end_LL   = squeeze(Uhis_LL(end, :, :));
Uz_LL      = U_end_LL(:, 3);

[uzmin_LL, idxLocal] = min(Uz_LL(bottomNodeIDs));
delta_LL             = -uzmin_LL;

Results.LL_only.name     = "LL_only (DC=0, PL=1.0)";
Results.LL_only.factors  = struct('DC', 0.0, 'PL', 1.0);
Results.LL_only.U_end    = U_end_LL;
Results.LL_only.Uhis     = Uhis_LL;
Results.LL_only.Uz       = Uz_LL;
Results.LL_only.delta    = delta_LL;
Results.LL_only.critNode = bottomNodeIDs(idxLocal);
Results.LL_only.uz_min   = uzmin_LL;
Results.LL_only.total_F  = -sum(nr_LL.load(:,4));

truss_strain_LL        = bar.Solve_Strain(node, U_end_LL);
Pu_LL                  = truss_strain_LL .* bar.E_vec .* bar.A_vec;
sigma_LL               = Pu_LL ./ bar.A_vec(:);
Results.LL_only.Pu     = Pu_LL;
Results.LL_only.sigma  = sigma_LL;

DO_PLOT_LL = false;
if DO_PLOT_LL
    plots.Plot_Shape_Bar_Stress(Results.LL_only.sigma);
end

fprintf("LL-only deflection (deck) = %.6f m at node %d  (span/%.0f)\n", ...
    Results.LL_only.delta, Results.LL_only.critNode, span_total / Results.LL_only.delta);

% Deflection limits (placeholder — verify per AASHTO 2.5.2.6.2)
delta_limit_360  = span_total / 360;
delta_limit_1000 = span_total / 1000;

fprintf("\n===== DEFLECTION CHECK (TODO: verify against AASHTO LRFD 2.5.2.6.2) =====\n");
fprintf("  span_total                        = %.3f m\n", span_total);
fprintf("  L_unit (module length, not limit) = %.3f m\n", L_unit);
fprintf("  Deflection under Service I        = %.6f m\n", delta_service);
fprintf("  Deflection under LL-only          = %.6f m\n", delta_LL);
fprintf("  Placeholder limit L/360           = %.6f m\n", delta_limit_360);
fprintf("  Placeholder limit L/1000          = %.6f m\n", delta_limit_1000);

if delta_LL <= delta_limit_360
    fprintf("  [PASS*] LL-only <= L/360 (placeholder): LL = L/%.0f\n", span_total/delta_LL);
else
    fprintf("  [FAIL*] LL-only >  L/360 (placeholder): LL = L/%.0f\n", span_total/delta_LL);
end

if delta_LL <= delta_limit_1000
    fprintf("  [PASS*] LL-only <= L/1000 (placeholder): LL = L/%.0f\n", span_total/delta_LL);
else
    fprintf("  [WARN*] LL-only >  L/1000 (placeholder): LL = L/%.0f\n", span_total/delta_LL);
end

% =========================================================================
% TODO (AASHTO Guide Spec for Pedestrian Bridges, Section 6):
%   Vibration/dynamic serviceability check not implemented.
%   Requires eigenvalue analysis: [V,D] = eig(K_global, M_global)
%   Vertical natural frequency must satisfy f_n >= 3.0 Hz (Guide Spec Sec.6)
%   or detailed dynamic analysis must be performed.
% =========================================================================

% Capacity-to-weight ratio
if exist('Cap','var') && Cap.converged
    Cap_to_W = Cap.total_F / Wbar;
else
    Cap_to_W = NaN;
end

fprintf("\n===== SYSTEM PERFORMANCE SUMMARY =====\n");
fprintf("  Bar self-weight (DC, bars only)  = %10.1f N\n",      Wbar);
fprintf("  Service I stiffness K            = %10.2e N/m\n",   K_service);
fprintf("  Stiffness-to-weight ratio        = %10.4f (N/m)/N\n", K_to_W);
if exist('Cap','var') && Cap.converged
    fprintf("  LRFD capacity (governing)        = %10.1f N\n",  Cap.total_F);
    fprintf("  Capacity-to-weight ratio         = %10.3f\n",    Cap_to_W);
else
    fprintf("  LRFD capacity                    = not computed\n");
end
fprintf("=======================================\n");

toc