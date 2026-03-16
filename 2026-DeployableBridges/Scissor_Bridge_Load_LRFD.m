clear all
close all
clc
tic

% =========================================================================
%  SCISSOR PEDESTRIAN BRIDGE - AASHTO LRFD STRUCTURAL ANALYSIS
%
%  Description : Full LRFD code check and capacity scan for a Scissor-
%                inspired pedestrian truss bridge under dead load (DC)
%                and pedestrian live load (PL).
%                The bridge includes passive truss bars (bar), standard
%                actuator bars (actBar), and 3-node / 4-node rotational
%                springs. All bar types are included in the DC self-weight.
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
%             Identical to Kirigami, Origami and Rolling bridges for
%             direct comparison.
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
% Verifies wall compactness; non-slender -> no effective area reduction (Q=1)
HSS = Check_HSS_Rect_Slenderness(bt, ht, barE, Fy);
fprintf("HSS local slenderness: b/t=%.2f, h/t=%.2f, lambda_r=%.2f -> %s\n", ...
        HSS.bt, HSS.ht, HSS.lambda_r, HSS.classStr);


%% -------------------------------------------------------------------------
%  SECTION 3: PANEL (SHEET) PROPERTIES
%             Panels are intentionally soft so that the truss carries all
%             global load; panel stiffness is set to ~1/1000 of steel.
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
%             Create empty containers for all element types.
%             Scissor bridge uses both bar (passive) and actBar (actuator),
%             plus 3-node and 4-node rotational springs.
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
%             - 8 nodes per section i = 1..N
%             - 4 closing nodes at x = L*N (end face)
%
%             Local layout per section (8 nodes):
%               1: (L*(i-1), 0, 0)   bottom-left start
%               2: (L*(i-1), L, 0)   bottom-right start
%               3: (L*(i-1), 0, L)   top-left start
%               4: (L*(i-1), L, L)   top-right start
%               5: (L*(i-1)+L/2, 0, L/2)   scissor apex-left
%               6: (L*(i-1)+L/2, L, L/2)   scissor apex-right
%               7: (L*(i-1)+L/2, 0, 0)     bottom-left mid
%               8: (L*(i-1)+L/2, L, 0)     bottom-right mid
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
        8*(i-1)+1, 8*(i-1)+2,  8*(i-1)+7;
        8*(i-1)+2, 8*(i-1)+7,  8*(i-1)+8;
        8*(i-1)+7, 8*(i-1)+8,  8*(i-1)+9;
        8*(i-1)+8, 8*(i-1)+10, 8*(i-1)+9;
        ];
end

cstNum    = size(cst.node_ijk_mat, 1);
cst.t_vec = panel_t * ones(cstNum, 1);
cst.E_vec = panel_E * ones(cstNum, 1);
cst.v_vec = panel_v * ones(cstNum, 1);

plots.Plot_Shape_CST_Number;


%% -------------------------------------------------------------------------
%  SECTION 8: BAR (TRUSS) ELEMENTS
%             Passive structural truss bars.
%             18 bars per section, plus 2 closing bars on the end face.
% -------------------------------------------------------------------------

for i = 1:N
    bar.node_ij_mat = [bar.node_ij_mat;
        % --- Bottom chord longitudinal bars (z = 0 plane) ---
        8*(i-1)+1, 8*(i-1)+7;      % bottom-left start to bottom-left mid
        8*(i-1)+7, 8*(i-1)+9;      % bottom-left mid to bottom-left end
        8*(i-1)+2, 8*(i-1)+8;      % bottom-right start to bottom-right mid
        8*(i-1)+8, 8*(i-1)+10;     % bottom-right mid to bottom-right end

        % --- Top chord and apex cross bars ---
        8*(i-1)+3, 8*(i-1)+4;      % top-left to top-right (end frame)
        8*(i-1)+9, 8*(i-1)+10;     % bottom-left end to bottom-right end (cross)
        8*(i-1)+1, 8*(i-1)+2;      % bottom-left start to bottom-right start (cross)
        8*(i-1)+7, 8*(i-1)+8;      % bottom-left mid to bottom-right mid (cross)

        % --- Scissor diagonal bars (left side) ---
        8*(i-1)+1, 8*(i-1)+5;      % bottom-left start to apex-left
        8*(i-1)+3, 8*(i-1)+5;      % top-left start to apex-left
        8*(i-1)+2, 8*(i-1)+6;      % bottom-right start to apex-right
        8*(i-1)+4, 8*(i-1)+6;      % top-right start to apex-right

        % --- Scissor diagonal bars (right side) ---
        8*(i-1)+5, 8*(i-1)+9;      % apex-left to bottom-left end
        8*(i-1)+5, 8*(i-1)+11;     % apex-left to top-left end
        8*(i-1)+6, 8*(i-1)+10;     % apex-right to bottom-right end
        8*(i-1)+6, 8*(i-1)+12;     % apex-right to top-right end

        % --- In-plane diagonal bars ---
        8*(i-1)+2, 8*(i-1)+7;      % bottom-right start to bottom-left mid
        8*(i-1)+8, 8*(i-1)+9;      % bottom-right mid to bottom-left end
        ];
end

% Closing bars on the terminal end face (x = L*N)
i = N + 1;
bar.node_ij_mat = [bar.node_ij_mat;
    8*(i-1)+1, 8*(i-1)+2;     % bottom-left to bottom-right, end face
    8*(i-1)+3, 8*(i-1)+4;     % top-left to top-right, end face
    ];

% Assign uniform cross-section properties to all passive bars
barNum_passive = size(bar.node_ij_mat, 1);
bar.A_vec      = barA * ones(barNum_passive, 1);
bar.E_vec      = barE * ones(barNum_passive, 1);

plots.Plot_Shape_Node_Number();
plots.Plot_Shape_Bar_Number();


%% -------------------------------------------------------------------------
%  SECTION 9: 3-NODE ROTATIONAL SPRINGS
%             Model scissor fold-line bending stiffness.
%             Stiffness = 100 * kspr (stiffened fold lines).
% -------------------------------------------------------------------------

for i = 1:N
    rotSpr3N.node_ijk_mat = [rotSpr3N.node_ijk_mat;
        % --- Left-side scissor fold springs ---
        8*(i-1)+1, 8*(i-1)+5,  8*(i-1)+11;
        8*(i-1)+3, 8*(i-1)+5,  8*(i-1)+9;
        8*(i-1)+2, 8*(i-1)+6,  8*(i-1)+12;
        8*(i-1)+4, 8*(i-1)+6,  8*(i-1)+10;

        % --- In-plane fold springs ---
        8*(i-1)+4, 8*(i-1)+3,  8*(i-1)+5;
        8*(i-1)+3, 8*(i-1)+4,  8*(i-1)+6;
        8*(i-1)+2, 8*(i-1)+1,  8*(i-1)+5;
        8*(i-1)+6, 8*(i-1)+2,  8*(i-1)+1;

        % --- Right-side scissor fold springs ---
        8*(i-1)+5, 8*(i-1)+11, 8*(i-1)+12;
        8*(i-1)+11, 8*(i-1)+12, 8*(i-1)+6;
        8*(i-1)+5, 8*(i-1)+9,  8*(i-1)+10;
        8*(i-1)+9, 8*(i-1)+10, 8*(i-1)+6;
        ];
end

rotNum3N = size(rotSpr3N.node_ijk_mat, 1);
rotSpr3N.rot_spr_K_vec = kspr * 100 * ones(rotNum3N, 1);

plots.Plot_Shape_Node_Number;
plots.Plot_Shape_RotSpr_3N_Number;


%% -------------------------------------------------------------------------
%  SECTION 10: 4-NODE ROTATIONAL SPRINGS
%              Applied at the bottom face panels for slab stability.
%              Stiffness 1e5 N·m/rad.
% -------------------------------------------------------------------------

for i = 1:N
    rotSpr4N.node_ijkl_mat = [rotSpr4N.node_ijkl_mat;
        8*(i-1)+1, 8*(i-1)+2, 8*(i-1)+7,  8*(i-1)+8;
        8*(i-1)+7, 8*(i-1)+8, 8*(i-1)+9,  8*(i-1)+10;
        ];
end

rotNum4N = size(rotSpr4N.node_ijkl_mat, 1);
rotSpr4N.rot_spr_K_vec = 1e5 * ones(rotNum4N, 1);

plots.Plot_Shape_RotSpr_4N_Number;


%% -------------------------------------------------------------------------
%  SECTION 11: ACTUATOR BAR ELEMENTS
%              Two batches of actBar:
%              Batch 1: vertical bars connecting bottom to top at each
%                       section face (N+1 faces x 2 bars = 2*(N+1) bars)
%              Batch 2: diagonal bars connecting bottom mid to apex at
%                       each section mid-plane (N x 2 bars)
%              All actBars use the same cross-section as passive bars.
% -------------------------------------------------------------------------

% Batch 1: vertical bars at section faces (including end face)
for i = 1:N
    actBar.node_ij_mat = [actBar.node_ij_mat;
        8*(i-1)+1, 8*(i-1)+3;    % bottom-left to top-left at section start
        8*(i-1)+2, 8*(i-1)+4;    % bottom-right to top-right at section start
        ];
end

% End-face vertical bars
i = N + 1;
actBar.node_ij_mat = [actBar.node_ij_mat;
    8*(i-1)+1, 8*(i-1)+3;    % bottom-left to top-left at end face
    8*(i-1)+2, 8*(i-1)+4;    % bottom-right to top-right at end face
    ];

% Batch 2: diagonal bars connecting mid-plane nodes to apex
for i = 1:N
    actBar.node_ij_mat = [actBar.node_ij_mat;
        8*(i-1)+7, 8*(i-1)+5;    % bottom-left mid to apex-left
        8*(i-1)+8, 8*(i-1)+6;    % bottom-right mid to apex-right
        ];
end

actBarNum     = size(actBar.node_ij_mat, 1);
actBar.A_vec  = barA * ones(actBarNum, 1);
actBar.E_vec  = barE * ones(actBarNum, 1);

plots.Plot_Shape_ActBar_Number;


%% -------------------------------------------------------------------------
%  SECTION 12: ASSEMBLY INITIALISATION  (build global stiffness structure)
% -------------------------------------------------------------------------

assembly.Initialize_Assembly;


%% -------------------------------------------------------------------------
%  SECTION 13: MATERIAL AND ANALYSIS CONTROL FLAGS
% -------------------------------------------------------------------------

rho_steel = 7850;   % Steel density (kg/m^3)
g         = 9.81;   % Gravitational acceleration (m/s^2)

DO_CODE_CHECK = true;   % PATH 1: run all load combinations + member DCR check
DO_CAPACITY   = true;   % PATH 2: scan load multiplier alpha to find capacity


%% =========================================================================
%  LOAD DEFINITION
%  AASHTO LRFD Section 3: Loads and Load Factors
% =========================================================================

%% -------------------------------------------------------------------------
%  STEP 1A: Identify deck (bottom chord) nodes for live load application
%           Bottom nodes are identified automatically as those at z = zMin
%           (z = 0). Consistent with Kirigami, Origami and Rolling bridges.
% -------------------------------------------------------------------------

Z_coords      = node.coordinates_mat(:, 3);
zMin_deck     = min(Z_coords);
zTol_deck     = max(1e-9, 1e-6 * max(1, abs(zMin_deck)));
bottomNodeIDs = find(abs(Z_coords - zMin_deck) <= zTol_deck);


%% -------------------------------------------------------------------------
%  STEP 1B: Build individual load case vectors
%
%  DC - Dead load: Component and attachment self-weight.
%       Includes BOTH passive bar and actuator bar self-weight.
%       (AASHTO Table 3.4.1-2: gamma_p_max = 1.25, gamma_p_min = 0.90)
%
%  PL - Pedestrian live load
%       AASHTO Section 3.6.1.6: w = 3.6 kPa for pedestrian bridges.
%       NOTE: qPL is currently set to 4.3 kPa (placeholder) and should
%       be updated to 3.6 kPa for full AASHTO compliance.
% -------------------------------------------------------------------------

LoadCase = struct();

nodeNum = size(node.coordinates_mat, 1);

% DC: passive bar self-weight
LoadCase_DC_passive = Build_LoadCase_DC_Bars(bar.node_ij_mat, ...
                          node.coordinates_mat, barA, rho_steel, g);

% DC: actuator bar self-weight (same cross-section as passive bars)
LoadCase_DC_actbar  = Build_LoadCase_DC_Bars(actBar.node_ij_mat, ...
                          node.coordinates_mat, barA, rho_steel, g);

% Combine passive + actuator DC loads into a single nodeNum x 4 matrix
DC_combined        = zeros(nodeNum, 4);
DC_combined(:, 1)  = (1:nodeNum)';

for row = 1:size(LoadCase_DC_passive, 1)
    nid = LoadCase_DC_passive(row, 1);
    DC_combined(nid, 2:4) = DC_combined(nid, 2:4) + LoadCase_DC_passive(row, 2:4);
end
for row = 1:size(LoadCase_DC_actbar, 1)
    nid = LoadCase_DC_actbar(row, 1);
    DC_combined(nid, 2:4) = DC_combined(nid, 2:4) + LoadCase_DC_actbar(row, 2:4);
end

LoadCase.DC = DC_combined;

% PL: uniform pedestrian pressure converted to equivalent nodal forces
span  = L * N;      % Total bridge span (m)
width = L;          % Bridge deck width (m)
qPL   = 4.3e3;      % Pedestrian load intensity (N/m^2) -- placeholder, update to 3.6e3
LoadCase.PL = Build_LoadCase_PL_Uniform(bottomNodeIDs, qPL, span, width);

% Sanity check: print total loads
Wbar_check = -sum(LoadCase.DC(:, 4));
PL_check   = -sum(LoadCase.PL(:, 4));
fprintf('DC total (bar + actBar self-weight)  = %.2f N\n', Wbar_check);
fprintf('PL total (alpha = 1.0)               = %.2f N\n', PL_check);


%% =========================================================================
%  STEP 2: LOAD COMBINATIONS  (AASHTO LRFD Table 3.4.1-1)
%
%  Service I   : Normal operational loads, used for deflection check.
%                Factors: 1.0 DC + 1.0 PL
%
%  Strength Ia : Maximum DC case, compression members typically govern.
%                Factors: 1.25 DC + 1.75 PL
%
%  Strength Ib : Minimum DC case, checks if reduced self-weight combined
%                with full live load produces worse tension in some members.
%                Factors: 0.90 DC + 1.75 PL
% =========================================================================

Combos = struct();

Combos.Service_1.name    = "Service_1";
Combos.Service_1.factors = struct('DC', 1.00, 'PL', 1.00);

Combos.Strength_1a.name    = "Strength_1a";
Combos.Strength_1a.factors = struct('DC', 1.25, 'PL', 1.75);   % DC maximum

Combos.Strength_1b.name    = "Strength_1b";
Combos.Strength_1b.factors = struct('DC', 0.90, 'PL', 1.75);   % DC minimum


%% =========================================================================
%  STEP 3: SOLVER AND BOUNDARY CONDITIONS
%
%  Asymmetric support conditions:
%    Start face (x=0) : nodes 1 and 2 fully pinned (Ux=Uy=Uz=0)
%    End face (x=L*N) : nodes 8N+1 and 8N+2 roller (Uy=Uz=0, Ux free)
%  This allows thermal expansion along the bridge axis (x-direction)
%  and is the physically correct condition for a simply supported bridge.
% =========================================================================

nr          = Solver_NR_Loading;
nr.assembly = assembly;

% All nodes initialised as free; then define supports
nodeNumVec = (1:nodeNum)';
nr.supp    = [nodeNumVec, zeros(nodeNum, 1), zeros(nodeNum, 1), zeros(nodeNum, 1)];

nr.supp(1,     2:4) = [1, 1, 1];   % start face, bottom-left: fully pinned
nr.supp(2,     2:4) = [1, 1, 1];   % start face, bottom-right: fully pinned
nr.supp(8*N+1, 2:4) = [0, 1, 1];   % end face, bottom-left: roller (Ux free)
nr.supp(8*N+2, 2:4) = [0, 1, 1];   % end face, bottom-right: roller (Ux free)


%% =========================================================================
%  STEP 4: MEMBER BUCKLING PARAMETERS  (shared by PATH 1 and PATH 2)
%
%  Effective length factor K = 1.0 for all members (pin-pin assumption).
%  Weak-axis radius of gyration ry governs (ry < rx for HSS 4x3).
%  KL/r verified to be well below the AISC limit of 200.
%
%  Note: Only passive bar members are checked in PATH 1 and PATH 2.
%  Actuator bars (actBar) should ideally be checked separately.
% =========================================================================

barNum_global  = size(bar.node_ij_mat, 1);

% Effective length factor (K = 1.0: both ends pinned)
Keff_vec       = ones(barNum_global, 1);

% Undeformed bar lengths from assembly
L0_vec_global  = bar.L0_vec(:);
KL_vec_global  = Keff_vec .* L0_vec_global;   % Effective lengths KL (m)

% Weak-axis radius of gyration applied uniformly to all members
r_vec_global   = ry * ones(barNum_global, 1);  % ry = 0.0287 m


%% =========================================================================
%  STEP 5: MEMBER TYPE CLASSIFICATION
%          Members are classified by the z-coordinates of their end nodes:
%            TopChord    : both nodes at z = zMax
%            BottomChord : both nodes at z = zMin
%            Web         : all other members
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
%          For each load combination:
%            1. Assemble factored nodal load vector
%            2. Solve for displacements (Newton-Raphson)
%            3. Recover bar forces (passive bars only)
%            4. Compute LRFD demand-to-capacity ratio (DCR = |Pu| / phiPn)
%               for tension yielding and compression buckling
%            5. Report maximum DCR and governing member
% =========================================================================

Results = struct();

if DO_CODE_CHECK

    comboNames = fieldnames(Combos);

    for c = 1:numel(comboNames)

        combo = Combos.(comboNames{c});
        fprintf("\n=== Running combination: %s ===\n", combo.name);

        % Solver settings
        nr.increStep = 1;
        nr.iterMax   = 50;
        nr.tol       = 1e-5;

        % Assemble factored load vector for this combination
        nr.load  = Build_Combo_Load(LoadCase, combo.factors, nodeNum);
        total_F  = -sum(nr.load(:, 4));
        fprintf("  Total factored vertical load = %.2f N\n", total_F);

        % Solve for nodal displacements
        Uhis    = nr.Solve;
        U_end   = squeeze(Uhis(end, :, :));    % Final displacement [nodeNum x 3]
        Uz      = U_end(:, 3);
        maxDisp = -min(Uz);                    % Maximum downward displacement (m)
        fprintf("  Max downward displacement    = %.6f m\n", maxDisp);

        % Recover axial forces in all passive bar members
        truss_strain   = bar.Solve_Strain(node, U_end);
        internal_force = truss_strain .* bar.E_vec .* bar.A_vec;   % Pu (N)

        barNum  = numel(internal_force);
        Pn      = NaN(barNum, 1);
        phi     = NaN(barNum, 1);
        phiPn   = NaN(barNum, 1);
        DCR     = NaN(barNum, 1);
        modeStr = cell(barNum, 1);
        passYN  = false(barNum, 1);

        % LRFD member check: tension yielding or compression buckling
        % phi*Pn per AISC 360 Chapter D (tension) and Chapter E (compression)
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

        % Guard against degenerate results
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

        % Report critical member geometry and forces
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

        % Store all results for this combination
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
%
%  The factored live load is scaled by a multiplier alpha:
%    DC_max case:  1.25*DC + alpha*(1.75*PL)
%    DC_min case:  0.90*DC + alpha*(1.75*PL)
%
%  Alpha is increased in coarse steps (0.5) until DCR >= 1.0 is detected,
%  then refined in fine steps (0.05) to locate the precise capacity point.
%  The case with the lower total capacity load governs design.
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

    Cap_all = struct();   % Stores capacity results for both DC cases

    for dc_idx = 1:numel(DC_cases)

        DC_factor = DC_cases(dc_idx);
        DC_label  = DC_labels(dc_idx);

        fprintf("\n--- Scanning %s ---\n", DC_label);

        % Initialise capacity result struct for this DC case
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

        % --- (1) Coarse scan: increment alpha until DCR first exceeds 1.0 ---
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

        % --- (2) Fine refinement within the bracketed interval ---
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

                % Print progress for each refinement step
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

        % Store result under appropriate field name
        if DC_factor == 1.25
            Cap_all.DC_max = Cap;
        else
            Cap_all.DC_min = Cap;
        end

    end   % end DC case loop

    % --- (3) Compare both DC cases and identify governing capacity ---
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

    % Select governing case: lower total factored load at DCR = 1.0 governs
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
%  Select which load combination to visualise.
%  Available options: "Service_1", "Strength_1a", "Strength_1b"
% =========================================================================

PLOT_CASE = "Strength_1a";

if isfield(Results, PLOT_CASE)

    DCR_plot   = Results.(PLOT_CASE).DCR;
    pass_plot  = Results.(PLOT_CASE).passYN;
    Pu_plot    = Results.(PLOT_CASE).Pu;
    A_plot     = bar.A_vec(:);
    sigma_plot = Pu_plot ./ A_plot;    % Axial stress (Pa); positive = tension

    [maxDCR_plot, kcrit_plot] = max(DCR_plot, [], 'omitnan');
    KLr_crit = KL_vec_global(kcrit_plot) / r_vec_global(kcrit_plot);

    fprintf("\n=== PLOTTING CASE: %s ===\n", PLOT_CASE);
    fprintf("  Max D/C = %.3f at bar #%d (%s)\n", ...
            maxDCR_plot, kcrit_plot, Results.(PLOT_CASE).critType);
    fprintf("  KL/r of governing member = %.2f  (limit = 200)\n", KLr_crit);

    % Uncomment to plot demand-to-capacity ratios across all members:
    % plots.Plot_Shape_Bar_DCR(DCR_plot);

    plots.Plot_Shape_Bar_Stress(sigma_plot);    % Axial stress distribution
    plots.Plot_Shape_Bar_Failure(pass_plot);    % Pass/fail map

else
    fprintf("\n[WARN] Results struct does not contain combination '%s'. Skipping plots.\n", PLOT_CASE);
end


%% =========================================================================
%  SYSTEM-LEVEL PERFORMANCE METRICS
%
%  Reports:
%    (1) Service I deflection check  (AASHTO Section 2.5.2.6.2)
%    (2) System stiffness and weight efficiency metrics
%    (3) Structural capacity relative to self-weight
% =========================================================================

Wbar = Wbar_check;   % Total bar self-weight (passive + actuator bars) (N)

% Service I displacement and stiffness
delta_service = Results.Service_1.maxDisp;
F_service     = Results.Service_1.total_F;
K_service     = F_service / delta_service;    % System stiffness (N/m)
K_to_W        = K_service / Wbar;

% -------------------------------------------------------------------------
%  DEFLECTION CHECK  (AASHTO Section 2.5.2.6.2)
%
%  Live load deflection limit for pedestrian bridges:
%    Standard case       : delta_LL <= L/360
%    Vibration-sensitive : delta_LL <= L/1000
%
%  Note: delta_service here conservatively includes DC contribution.
%  For a strict LL-only check, solve with factors = {DC:0, PL:1.0}.
% -------------------------------------------------------------------------

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

% Capacity-to-weight ratio
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