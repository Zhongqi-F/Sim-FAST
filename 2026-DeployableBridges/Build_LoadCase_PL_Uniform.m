function loadPL = Build_LoadCase_PL_Uniform(bottomNodeIDs, qPL, span, width)
% qPL: pedestrian area load in N/m^2
% span: bridge length (m), width: deck width (m)

Adeck = span * width;
Ftot  = qPL * Adeck;                       % total vertical load (N)
fn    = -Ftot / numel(bottomNodeIDs);      % nodal vertical load (N), downward is negative

loadPL = [bottomNodeIDs(:), ...
          zeros(numel(bottomNodeIDs),1), ...
          zeros(numel(bottomNodeIDs),1), ...
          fn * ones(numel(bottomNodeIDs),1)];
end