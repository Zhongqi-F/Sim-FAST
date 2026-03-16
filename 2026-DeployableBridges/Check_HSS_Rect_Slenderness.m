function HSS = Check_HSS_Rect_Slenderness(bt, ht, E, Fy, varargin)
% Check local slenderness for rectangular HSS walls under axial compression
% per AASHTO LRFD Section 6.9.4.2.1 (Table 6.9.4.2.1-1).
%
% Inputs:
%   bt : b/t ratio for the wall in the "B" direction (dimensionless)
%   ht : h/t ratio for the wall in the "H" direction (dimensionless)
%   E  : modulus of elasticity (Pa)
%   Fy : specified minimum yield strength (Pa)
%
% Optional name-value:
%   'Forming' : 'cold-formed' (default) or 'hot-formed'
%       - cold-formed rectangular HSS walls: lambda_r = 1.28*sqrt(E/Fy)
%       - hot-formed rectangular HSS walls (and some box members): 1.40*sqrt(E/Fy)
%
% Output struct HSS:
%   .bt, .ht, .lambda_r
%   .isSlender_b, .isSlender_h, .isSlender
%   .classStr
%   .forming, .coeff

    % -------------------------
    % Parse optional arguments
    % -------------------------
    p = inputParser;
    p.addParameter('Forming', 'cold-formed', @(s)ischar(s) || isstring(s));
    p.parse(varargin{:});
    forming = lower(string(p.Results.Forming));

    % -------------------------
    % Guard
    % -------------------------
    if bt <= 0 || ht <= 0 || E <= 0 || Fy <= 0
        error("Invalid inputs: require bt>0, ht>0, E>0, Fy>0.");
    end

    % -------------------------
    % AASHTO Table 6.9.4.2.1-1 limits
    % -------------------------
    switch forming
        case "cold-formed"
            coeff = 1.28;   % Walls of square/rectangular cold-formed HSS
        case "hot-formed"
            coeff = 1.40;   % Walls of square/rectangular hot-formed HSS (and some box cases)
        otherwise
            error("Unknown Forming option. Use 'cold-formed' or 'hot-formed'.");
    end

    lambda_r = coeff * sqrt(E / Fy);

    % -------------------------
    % Classification
    % -------------------------
    HSS.bt       = bt;
    HSS.ht       = ht;
    HSS.lambda_r = lambda_r;
    HSS.forming  = forming;
    HSS.coeff    = coeff;

    HSS.isSlender_b = (bt > lambda_r);
    HSS.isSlender_h = (ht > lambda_r);
    HSS.isSlender   = (HSS.isSlender_b || HSS.isSlender_h);

    if HSS.isSlender
        HSS.classStr = "Slender per AASHTO 6.9.4.2.1 (local buckling effects may apply)";
    else
        HSS.classStr = "Non-slender per AASHTO 6.9.4.2.1 (local buckling effects neglected)";
    end
end