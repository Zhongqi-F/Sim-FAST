clear all
close all
clc
tic

% =========================================================================
%  SCISSOR PEDESTRIAN BRIDGE 2 - AASHTO LRFD STRUCTURAL ANALYSIS
%
%  Description : Full LRFD code check and capacity scan for a Scissor V2
%                pedestrian truss bridge under dead load (DC) and
%                pedestrian live load (PL).
%                V2 extends V1 by adding two top-face mid-plane nodes
%                (local 9, 10 at z=L) per section, enriching the top
%                chord connectivity and rotational spring coverage.
%                Passive bars (bar), standard actuator bars (actBar),
%                and 3-node / 4-node rotational springs are all included.
%                Supports are asymmetric: start face fully pinned,
%                end face roller (x-direction free).
%
%  Load cases  : Service I   (1.0  DC + 1.0  PL)
%                Strength Ia (1.25 DC + 1.75 PL)   <- DC maximum
%                Strength Ib (0.90 DC + 1.75 PL)   <- DC minimum
%
%  Code refs   : AASHTO LRFD Bridge Design Specifications
%                AISC 360 (member strength, compression buckling)
%
%  Units       : SI throughout (N, m, Pa)
% =========================================================================


%% -------------------------------------------------------------------------
%  SECTION 1: GEOMETRY PARAMETERS
% -------------------------------------------------------------------------

N = 8;      % Number of repeating sections along bridge span
H = 2;      % Bridge height (m)
L = 2;      % Section length = bridge width (m)


%% -------------------------------------------------------------------------
%  SECTION 2: CROSS-SECTION PROPERTIES  (HSS 4x3x5/16, A500 Grade C)
%             Source: AISC Steel Construction Manual
%             Identical to Kirigami, Origami, Rolling and Scissor V1
%             bridges for direct comparison.
% -------------------------------------------------------------------------

barE  = 200e9;      % Elastic modulus, steel (Pa)
Fy    = 345e6;      % Yield strength, A500 Gr C ~ Grade 50 (Pa)

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
HSS = Check_HSS_Rect_Slenderness(bt, ht, barE, Fy);
fprintf("HSS local slenderness: b/t=%.2f, h/t=%.2f, lambda_r=%.2f -> %s\n", ...
        HSS.bt, HSS.ht, HSS.lambda_r, HSS.classStr);


%% -------------------------------------------------------------------------
%  SECTION 3: PANEL (SHEET) PROPERTIES
% -------------------------------------------------------------------------

panel_E = 2e8;      % Panel elastic modulus (Pa)  ~200 MPa (soft)
panel_t = 0.01;     % Panel thickness (m)
panel_v = 0.3;      % Panel Poisson's ratio


%% -------------------------------------------------------------------------
%  SECTION 4: ROTATIONAL SPRING STIFFNESS
%             3-node springs model scissor fold-line bending.
%             Stiffness derived from EI/L of the diagonal bar.
% -------------------------------------------------------------------------

barL = sqrt(H^2 + L^2);         % Diagonal bar length (m)
I    = ry^2 * barA;              % Weak-axis moment of inertia (m^4)
kspr = barE * I / barL;          % Rotational spring stiffness (N·m/rad)


%% -------------------------------------------------------------------------
%  SECTION 5: ASSEMBLY INITIALISATION
% -------------------------------------------------------------------------

assembly  = Assembly_Scissor_Bridge();
node      = Elements_Nodes();
cst       = Vec_Elements_CST;
rotSpr3N  = CD_Elements_RotSprings_3N;
rotSpr4N  = Vec_Elements_RotSprings_4N;
bar       = Vec_Elements_Bars;
actBar    = Std_Elements_Bars;

assembly.cst        = cst;
assembly.node       = node;
assembly.bar        = bar;
assembly.rot_spr_3N = rotSpr3N;
assembly.rot_spr_4N = rotSpr4N;
assembly.actBar     = actBar;


%% -------------------------------------------------------------------------
%  SECTION 6: NODE COORDINATES
%             Global node numbering:
%             - 10 nodes per section i = 1..N
%             - 4 closing nodes at x = L*N (end face)
%
%             Local layout per section (10 nodes):
%               1: (L*(i-1), 0, 0)        bottom-left start
%               2: (L*(i-1), L, 0)        bottom-right start
%               3: (L*(i-1), 0, L)        top-left start
%               4: (L*(i-1), L, L)        top-right start
%               5: (L*(i-1)+L/2, 0, L/2)  scissor apex-left
%               6: (L*(i-1)+L/2, L, L/2)  scissor apex-right
%               7: (L*(i-1)+L/2, 0, 0)    bottom-left mid
%               8: (L*(i-1)+L/2, L, 0)    bottom-right mid
%               9: (L*(i-1)+L/2, 0, L)    top-left mid      [V2 addition]
%              10: (L*(i-1)+L/2, L, L)    top-right mid     [V2 addition]
% -------------------------------------------------------------------------

for i = 1:N
    node.coordinates_mat = [node.coordinates_mat;
        L*(i-1),     0, 0;      % local 1: bottom-left start
        L*(i-1),     L, 0;      % local 2: bottom-right start
        L*(i-1),     0, L;      % local 3: top-left start
        L*(i-1),     L, L;      % local 4: top-right start
        (i-1)*L+L/2, 0, L/2;   % local 5: scissor apex-left
        (i-1)*L+L/2, L, L/2;   % local 6: scissor apex-right
        (i-1)*L+L/2, 0, 0;     % local 7: bottom-left mid
        (i-1)*L+L/2, L, 0;     % local 8: bottom-right mid
        (i-1)*L+L/2, 0, L;     % local 9: top-left mid  [V2]
        (i-1)*L+L/2, L, L;     % local 10: top-right mid [V2]
        ];
end

% End-face closing nodes (x = L*N)
node.coordinates_mat = [node.coordinates_mat;
    L*N, 0, 0;      % closing node 1: bottom-left
    L*N, L, 0;      % closing node 2: bottom-right
    L*N, 0, L;      % closing node 3: top-left
    L*N, L, L;      % closing node 4: top-right
    ];


%% -------------------------------------------------------------------------
%  SECTION 7: PLOTTING SETUP
% -------------------------------------------------------------------------

plots = Plot_Scissor_Bridge;
plots.assembly    = assembly;
plots.displayRange = [-1; 2*N+1; -1; 3; -1; 3];
plots.viewAngle1  = 20;
plots.viewAngle2  = 20;

plots.Plot_Shape_Node_Number;


%% -------------------------------------------------------------------------
%  SECTION 8: CST PANEL ELEMENTS
%             4 triangular panels per section on the bottom face.
% -------------------------------------------------------------------------

for i = 1:N
    cst.node_ijk_mat = [cst.node_ijk_mat;
        10*(i-1)+1,  10*(i-1)+2,  10*(i-1)+7;
        10*(i-1)+2,  10*(i-1)+7,  10*(i-1)+8;
        10*(i-1)+7,  10*(i-1)+8,  10*(i-1)+11;
        10*(i-1)+8,  10*(i-1)+12, 10*(i-1)+11;
        ];
end

cstNum    = size(cst.node_ijk_mat, 1);
cst.t_vec = panel_t * ones(cstNum, 1);
cst.E_vec = panel_E * ones(cstNum, 1);
cst.v_vec = panel_v * ones(cstNum, 1);

plots.Plot_Shape_CST_Number;


%% -------------------------------------------------------------------------
%  SECTION 9: BAR (TRUSS) ELEMENTS
%             Passive structural truss bars.
%             24 bars per section (V1 had 18; V2 adds 6 top-chord bars
%             connecting local nodes 9, 10 to their neighbours).
%             Plus 2 closing bars on the end face.
% -------------------------------------------------------------------------

for i = 1:N
    bar.node_ij_mat = [bar.node_ij_mat;
        % --- Bottom chord longitudinal bars (z = 0 plane) ---
        10*(i-1)+1,  10*(i-1)+7;      % bottom-left start to bottom-left mid
        10*(i-1)+7,  10*(i-1)+11;     % bottom-left mid to bottom-left end
        10*(i-1)+2,  10*(i-1)+8;      % bottom-right start to bottom-right mid
        10*(i-1)+8,  10*(i-1)+12;     % bottom-right mid to bottom-right end

        % --- Top chord longitudinal bars (z = L plane) [V2 enriched] ---
        10*(i-1)+3,  10*(i-1)+9;      % top-left start to top-left mid
        10*(i-1)+9,  10*(i-1)+13;     % top-left mid to top-left end
        10*(i-1)+4,  10*(i-1)+10;     % top-right start to top-right mid
        10*(i-1)+10, 10*(i-1)+14;     % top-right mid to top-right end

        % --- Top and bottom cross bars ---
        10*(i-1)+3,  10*(i-1)+4;      % top-left to top-right (start face)
        10*(i-1)+9,  10*(i-1)+10;     % top-left mid to top-right mid [V2]
        10*(i-1)+1,  10*(i-1)+2;      % bottom-left to bottom-right (start face)
        10*(i-1)+7,  10*(i-1)+8;      % bottom-left mid to bottom-right mid

        % --- Scissor diagonal bars (left side) ---
        10*(i-1)+1,  10*(i-1)+5;      % bottom-left start to apex-left
        10*(i-1)+3,  10*(i-1)+5;      % top-left start to apex-left
        10*(i-1)+2,  10*(i-1)+6;      % bottom-right start to apex-right
        10*(i-1)+4,  10*(i-1)+6;      % top-right start to apex-right

        % --- Scissor diagonal bars (right side) ---
        10*(i-1)+5,  10*(i-1)+13;     % apex-left to top-left end
        10*(i-1)+5,  10*(i-1)+11;     % apex-left to bottom-left end
        10*(i-1)+6,  10*(i-1)+12;     % apex-right to bottom-right end
        10*(i-1)+6,  10*(i-1)+14;     % apex-right to top-right end

        % --- Top chord diagonals [V2] ---
        10*(i-1)+3,  10*(i-1)+10;     % top-left start to top-right mid
        10*(i-1)+13, 10*(i-1)+10;     % top-left end to top-right mid

        % --- In-plane diagonal bars ---
        10*(i-1)+2,  10*(i-1)+7;      % bottom-right start to bottom-left mid
        10*(i-1)+8,  10*(i-1)+11;     % bottom-right mid to bottom-left end
        ];
end

% Closing bars on the terminal end face (x = L*N)
i = N + 1;
bar.node_ij_mat = [bar.node_ij_mat;
    10*(i-1)+1, 10*(i-1)+2;     % bottom-left to bottom-right, end face
    10*(i-1)+3, 10*(i-1)+4;     % top-left to top-right, end face
    ];

% Assign uniform cross-section properties to all passive bars
barNum_passive = size(bar.node_ij_mat, 1);
bar.A_vec      = barA * ones(barNum_passive, 1);
bar.E_vec      = barE * ones(barNum_passive, 1);

plots.Plot_Shape_Node_Number();
plots.Plot_Shape_Bar_Number();


%% -------------------------------------------------------------------------
%  SECTION 10: 3-NODE ROTATIONAL SPRINGS
%              16 springs per section (V1 had 12; V2 adds 4 for top mid).
%              Stiffness = 100 * kspr.
% -------------------------------------------------------------------------

for i = 1:N
    rotSpr3N.node_ijk_mat = [rotSpr3N.node_ijk_mat;
        % --- Scissor fold springs ---
        10*(i-1)+1,  10*(i-1)+5,  10*(i-1)+13;
        10*(i-1)+3,  10*(i-1)+5,  10*(i-1)+11;
        10*(i-1)+2,  10*(i-1)+6,  10*(i-1)+14;
        10*(i-1)+4,  10*(i-1)+6,  10*(i-1)+12;

        % --- In-plane fold springs ---
        10*(i-1)+4,  10*(i-1)+3,  10*(i-1)+5;
        10*(i-1)+3,  10*(i-1)+4,  10*(i-1)+6;
        10*(i-1)+2,  10*(i-1)+1,  10*(i-1)+5;
        10*(i-1)+6,  10*(i-1)+2,  10*(i-1)+1;

        % --- Bottom-face right-side fold springs ---
        10*(i-1)+5,  10*(i-1)+11, 10*(i-1)+12;
        10*(i-1)+11, 10*(i-1)+12, 10*(i-1)+6;

        % --- Top-face left-side fold springs [V2 enriched] ---
        10*(i-1)+5,  10*(i-1)+13, 10*(i-1)+14;
        10*(i-1)+13, 10*(i-1)+14, 10*(i-1)+6;

        % --- Top mid-plane crease springs [V2 addition] ---
        10*(i-1)+11, 10*(i-1)+13, 10*(i-1)+14;
        10*(i-1)+13, 10*(i-1)+14, 10*(i-1)+12;
        10*(i-1)+14, 10*(i-1)+12, 10*(i-1)+11;
        10*(i-1)+12, 10*(i-1)+11, 10*(i-1)+13;
        ];
end

rotNum3N = size(rotSpr3N.node_ijk_mat, 1);
rotSpr3N.rot_spr_K_vec = kspr * 100 * ones(rotNum3N, 1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;


%% -------------------------------------------------------------------------
%  SECTION 11: 4-NODE ROTATIONAL SPRINGS
%              4 springs per section (V1 had 2; V2 adds top-face springs).
%              Stiffness 1e5 N·m/rad.
% -------------------------------------------------------------------------

for i = 1:N
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
        % Bottom face springs
        10*(i-1)+1,  10*(i-1)+2,  10*(i-1)+7,  10*(i-1)+8;
        10*(i-1)+7,  10*(i-1)+8,  10*(i-1)+11, 10*(i-1)+12;
        % Top face springs [V2 addition]
        10*(i-1)+4,  10*(i-1)+3,  10*(i-1)+10, 10*(i-1)+9;
        10*(i-1)+10, 10*(i-1)+9,  10*(i-1)+14, 10*(i-1)+13;
        ];
end

rotNum4N = size(rotSpr4N.node_ijkl_mat, 1);
rotSpr4N.rot_spr_K_vec = 1e5 * ones(rotNum4N, 1);

plots.Plot_Shape_RotSpr_4N_Number;


%% -------------------------------------------------------------------------
%  SECTION 12: ACTUATOR BAR ELEMENTS
%              Batch 1: vertical bars at section faces (N+1 faces x 2)
%              Batch 2: diagonal bars at mid-plane [V2: 4 per section
%                       vs V1's 2; adds apex-to-top-mid connections]
% -------------------------------------------------------------------------

% Batch 1: vertical bars at section faces
for i = 1:N
    actBar.node_ij_mat = [actBar.node_ij_mat;
        10*(i-1)+1, 10*(i-1)+3;    % bottom-left to top-left at section start
        10*(i-1)+2, 10*(i-1)+4;    % bottom-right to top-right at section start
        ];
end

% End-face vertical bars
i = N + 1;
actBar.node_ij_mat = [actBar.node_ij_mat;
    10*(i-1)+1, 10*(i-1)+3;
    10*(i-1)+2, 10*(i-1)+4;
    ];

% Batch 2: mid-plane diagonal bars [V2: 4 per section]
for i = 1:N
    actBar.node_ij_mat = [actBar.node_ij_mat;
        10*(i-1)+7, 10*(i-1)+5;    % bottom-left mid to apex-left
        10*(i-1)+5, 10*(i-1)+9;    % apex-left to top-left mid  [V2]
        10*(i-1)+6, 10*(i-1)+10;   % apex-right to top-right mid [V2]
        10*(i-1)+8, 10*(i-1)+6;    % bottom-right mid to apex-right
        ];
end

actBarNum     = size(actBar.node_ij_mat, 1);
actBar.A_vec  = barA * ones(actBarNum, 1);
actBar.E_vec  = barE * ones(actBarNum, 1);

plots.Plot_Shape_ActBar_Number;


%% -------------------------------------------------------------------------
%  SECTION 13: ASSEMBLY INITIALISATION
% -------------------------------------------------------------------------

assembly.Initialize_Assembly;


%% -------------------------------------------------------------------------
%  SECTION 14: MATERIAL AND ANALYSIS CONTROL FLAGS
% -------------------------------------------------------------------------

rho_steel = 7850;   % Steel density (kg/m^3)
g         = 9.81;   % Gravitational acceleration (m/s^2)

DO_CODE_CHECK = true;
DO_CAPACITY   = true;


%% =========================================================================
%  LOAD DEFINITION
%  AASHTO LRFD Section 3: Loads and Load Factors
% =========================================================================

%% -------------------------------------------------------------------------
%  STEP 1A: Identify deck (bottom chord) nodes for live load application
%           z = zMin (z = 0) auto-detection, consistent with all other bridges.
% -------------------------------------------------------------------------

Z_coords      = node.coordinates_mat(:, 3);
zMin_deck     = min(Z_coords);
zTol_deck     = max(1e-9, 1e-6 * max(1, abs(zMin_deck)));
bottomNodeIDs = find(abs(Z_coords - zMin_deck) <= zTol_deck);


%% -------------------------------------------------------------------------
%  STEP 1B: Build individual load case vectors
%
%  DC - Dead load: passive bar + actuator bar self-weight combined.
%  PL - Pedestrian live load (4.3 kPa placeholder; update to 3.6 kPa).
% -------------------------------------------------------------------------

LoadCase = struct();
nodeNum  = size(node.coordinates_mat, 1);

% DC: passive bar self-weight
LoadCase_DC_passive = Build_LoadCase_DC_Bars(bar.node_ij_mat, ...
                          node.coordinates_mat, barA, rho_steel, g);

% DC: actuator bar self-weight
LoadCase_DC_actbar  = Build_LoadCase_DC_Bars(actBar.node_ij_mat, ...
                          node.coordinates_mat, barA, rho_steel, g);

% Combine into single nodeNum x 4 matrix
DC_combined       = zeros(nodeNum, 4);
DC_combined(:, 1) = (1:nodeNum)';

for row = 1:size(LoadCase_DC_passive, 1)
    nid = LoadCase_DC_passive(row, 1);
    DC_combined(nid, 2:4) = DC_combined(nid, 2:4) + LoadCase_DC_passive(row, 2:4);
end
for row = 1:size(LoadCase_DC_actbar, 1)
    nid = LoadCase_DC_actbar(row, 1);
    DC_combined(nid, 2:4) = DC_combined(nid, 2:4) + LoadCase_DC_actbar(row, 2:4);
end

LoadCase.DC = DC_combined;

% PL: uniform pedestrian pressure
span  = L * N;      % Total bridge span (m)
width = L;          % Bridge deck width (m)
qPL   = 4.3e3;      % Pedestrian load intensity (N/m^2) -- placeholder, update to 3.6e3
LoadCase.PL = Build_LoadCase_PL_Uniform(bottomNodeIDs, qPL, span, width);

Wbar_check = -sum(LoadCase.DC(:, 4));
PL_check   = -sum(LoadCase.PL(:, 4));
fprintf('DC total (bar + actBar self-weight)  = %.2f N\n', Wbar_check);
fprintf('PL total (alpha = 1.0)               = %.2f N\n', PL_check);


%% =========================================================================
%  STEP 2: LOAD COMBINATIONS  (AASHTO LRFD Table 3.4.1-1)
% =========================================================================

Combos = struct();

Combos.Service_1.name    = "Service_1";
Combos.Service_1.factors = struct('DC', 1.00, 'PL', 1.00);

Combos.Strength_1a.name    = "Strength_1a";
Combos.Strength_1a.factors = struct('DC', 1.25, 'PL', 1.75);

Combos.Strength_1b.name    = "Strength_1b";
Combos.Strength_1b.factors = struct('DC', 0.90, 'PL', 1.75);


%% =========================================================================
%  STEP 3: SOLVER AND BOUNDARY CONDITIONS
%
%  Asymmetric supports (same as Scissor V1):
%    Start face: nodes 1, 2 fully pinned (Ux=Uy=Uz=0)
%    End face  : nodes 10*N+1, 10*N+2 roller (Uy=Uz=0, Ux free)
% =========================================================================

nr          = Solver_NR_Loading;
nr.assembly = assembly;

nodeNumVec = (1:nodeNum)';
nr.supp    = [nodeNumVec, zeros(nodeNum,1), zeros(nodeNum,1), zeros(nodeNum,1)];

nr.supp(1,      2:4) = [1, 1, 1];   % start face, bottom-left: fully pinned
nr.supp(2,      2:4) = [1, 1, 1];   % start face, bottom-right: fully pinned
nr.supp(10*N+1, 2:4) = [0, 1, 1];   % end face, bottom-left: roller (Ux free)
nr.supp(10*N+2, 2:4) = [0, 1, 1];   % end face, bottom-right: roller (Ux free)


%% =========================================================================
%  STEP 4: MEMBER BUCKLING PARAMETERS
% =========================================================================

barNum_global  = size(bar.node_ij_mat, 1);
Keff_vec       = ones(barNum_global, 1);
L0_vec_global  = bar.L0_vec(:);
KL_vec_global  = Keff_vec .* L0_vec_global;
r_vec_global   = ry * ones(barNum_global, 1);


%% =========================================================================
%  STEP 5: MEMBER TYPE CLASSIFICATION
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
    sum(memberType_global == "TopChord"),    ...
    sum(memberType_global == "BottomChord"), ...
    sum(memberType_global == "Web"));


%% =========================================================================
%  PATH 1: PER-COMBINATION CODE CHECK
% =========================================================================

Results = struct();

if DO_CODE_CHECK

    comboNames = fieldnames(Combos);

    for c = 1:numel(comboNames)

        combo = Combos.(comboNames{c});
        fprintf("\n=== Running combination: %s ===\n", combo.name);

        nr.increStep = 1;
        nr.iterMax   = 50;
        nr.tol       = 1e-5;

        nr.load  = Build_Combo_Load(LoadCase, combo.factors, nodeNum);
        total_F  = -sum(nr.load(:, 4));
        fprintf("  Total factored vertical load = %.2f N\n", total_F);

        Uhis    = nr.Solve;
        U_end   = squeeze(Uhis(end, :, :));
        maxDisp = -min(U_end(:, 3));
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
            [passi, modeStri, Pni, phii, phiPni, DCRi] = Check_Truss_LRFD( ...
                internal_force(k), bar.A_vec(k), bar.E_vec(k), ...
                KL_vec_global(k), r_vec_global(k), Fy);
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
            Results.(comboNames{c}).total_F = total_F;
            Results.(comboNames{c}).maxDisp = maxDisp;
            Results.(comboNames{c}).maxDCR  = maxDCR;
            Results.(comboNames{c}).passYN  = passYN;
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

        Results.(comboNames{c}).total_F   = total_F;
        Results.(comboNames{c}).maxDisp   = maxDisp;
        Results.(comboNames{c}).maxDCR    = maxDCR;
        Results.(comboNames{c}).critBar   = kcrit;
        Results.(comboNames{c}).critMode  = modeStr{kcrit};
        Results.(comboNames{c}).critType  = critType;
        Results.(comboNames{c}).critNodes = [n1, n2];
        Results.(comboNames{c}).critLen   = Lcrit;
        Results.(comboNames{c}).passYN    = passYN;
        Results.(comboNames{c}).DCR       = DCR;
        Results.(comboNames{c}).Pu        = internal_force;
        Results.(comboNames{c}).Pn        = Pn;
        Results.(comboNames{c}).phi       = phi;
        Results.(comboNames{c}).phiPn     = phiPn;
        Results.(comboNames{c}).modeStr   = modeStr;

    end
end


%% =========================================================================
%  PATH 2: LRFD CAPACITY SCAN  (Strength I, DC maximum and DC minimum)
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
                    KL_vec_global(k), r_vec_global(k), Fy);
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
                        KL_vec_global(k), r_vec_global(k), Fy);
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
    KLr_crit = KL_vec_global(kcrit_plot) / r_vec_global(kcrit_plot);

    fprintf("\n=== PLOTTING CASE: %s ===\n", PLOT_CASE);
    fprintf("  Max D/C = %.3f at bar #%d (%s)\n", ...
            maxDCR_plot, kcrit_plot, Results.(PLOT_CASE).critType);
    fprintf("  KL/r of governing member = %.2f  (limit = 200)\n", KLr_crit);

    plots.Plot_Shape_Bar_Stress(sigma_plot);
    plots.Plot_Shape_Bar_Failure(pass_plot);

else
    fprintf("\n[WARN] Results struct does not contain combination '%s'. Skipping plots.\n", PLOT_CASE);
end


%% =========================================================================
%  SYSTEM-LEVEL PERFORMANCE METRICS
% =========================================================================

Wbar = Wbar_check;

delta_service = Results.Service_1.maxDisp;
F_service     = Results.Service_1.total_F;
K_service     = F_service / delta_service;
K_to_W        = K_service / Wbar;

delta_limit_360  = span / 360;
delta_limit_1000 = span / 1000;

fprintf("\n===== SERVICE I DEFLECTION CHECK (AASHTO 2.5.2.6.2) =====\n");
fprintf("  Bridge span                   = %.3f m\n", span);
fprintf("  Max deflection (Service I)    = %.6f m  (L/%.0f)\n", ...
        delta_service, span / delta_service);
fprintf("  Limit L/360                   = %.6f m\n", delta_limit_360);
fprintf("  Limit L/1000 (vibration)      = %.6f m\n", delta_limit_1000);

if delta_service <= delta_limit_360
    fprintf("  [PASS] L/360  check: delta = L/%.0f\n", span / delta_service);
else
    fprintf("  [FAIL] L/360  check: delta = L/%.0f  exceeds L/360 limit!\n", ...
            span / delta_service);
end

if delta_service <= delta_limit_1000
    fprintf("  [PASS] L/1000 check: delta = L/%.0f\n", span / delta_service);
else
    fprintf("  [WARN] L/1000 check: delta = L/%.0f  exceeds L/1000 limit\n", ...
            span / delta_service);
end

if Cap.converged
    Cap_to_W = Cap.total_F / Wbar;
else
    Cap_to_W = NaN;
end

fprintf("\n===== SYSTEM PERFORMANCE SUMMARY =====\n");
fprintf("  Bar self-weight (passive+act) = %10.1f N\n",     Wbar);
fprintf("  Service I stiffness K         = %10.2e N/m\n",  K_service);
fprintf("  Stiffness-to-weight ratio     = %10.4f (N/m)/N\n", K_to_W);
fprintf("  LRFD capacity (governing)     = %10.1f N\n",     Cap.total_F);
fprintf("  Capacity-to-weight ratio      = %10.3f\n",       Cap_to_W);
fprintf("=======================================\n");

toc