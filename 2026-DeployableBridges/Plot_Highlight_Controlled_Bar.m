function h = Plot_Highlight_Controlled_Bar(nodeCoord, barIJ, barID, varargin)
%Plot_Highlight_Controlled_Bar  Plot truss geometry and highlight one bar.
%
% Usage:
%   Plot_Highlight_Controlled_Bar(node.coordinates_mat, bar.node_ij_mat, 144);
%
% Optional name-value pairs:
%   'ShowAllBars'  (true/false) default true
%   'AllColor'     [r g b]      default [0.75 0.75 0.75]
%   'AllLineWidth' scalar       default 0.8
%   'HiColor'      [r g b]      default [1 0 0]
%   'HiLineWidth'  scalar       default 3.0
%   'Label'        true/false   default true
%   'View'         [az el]      default [] (no change)
%   'Title'        char/string  default sprintf('Highlight bar #%d',barID)
%
% Output:
%   h struct with handles: h.fig, h.ax, h.allBars (optional), h.hiBar, h.txt

    p = inputParser;
    p.addParameter('ShowAllBars', true, @(x)islogical(x) || isnumeric(x));
    p.addParameter('AllColor', [0.75 0.75 0.75], @(x)isnumeric(x) && numel(x)==3);
    p.addParameter('AllLineWidth', 0.8, @(x)isnumeric(x) && isscalar(x));
    p.addParameter('HiColor', [1 0 0], @(x)isnumeric(x) && numel(x)==3);
    p.addParameter('HiLineWidth', 3.0, @(x)isnumeric(x) && isscalar(x));
    p.addParameter('Label', true, @(x)islogical(x) || isnumeric(x));
    p.addParameter('View', [], @(x)isnumeric(x) && (isempty(x) || numel(x)==2));
    p.addParameter('Title', "", @(x)ischar(x) || isstring(x));
    p.parse(varargin{:});

    showAll = logical(p.Results.ShowAllBars);
    allC    = p.Results.AllColor;
    allLW   = p.Results.AllLineWidth;
    hiC     = p.Results.HiColor;
    hiLW    = p.Results.HiLineWidth;
    doLabel = logical(p.Results.Label);
    vw      = p.Results.View;

    barNum = size(barIJ, 1);
    if barID < 1 || barID > barNum
        error("barID=%d out of range (1..%d).", barID, barNum);
    end

    fig = figure;
    ax = axes(fig);
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');

    ttl = p.Results.Title;
    if strlength(string(ttl))==0
        ttl = sprintf("Highlight bar #%d", barID);
    end
    title(ax, ttl);

    h = struct();
    h.fig = fig;
    h.ax  = ax;
    h.allBars = gobjects(0);
    h.hiBar   = gobjects(1);
    h.txt     = gobjects(1);

    % Plot all members (optional)
    if showAll
        % faster plotting: loop once
        h.allBars = gobjects(barNum,1);
        for e = 1:barNum
            n1 = barIJ(e,1); n2 = barIJ(e,2);
            p1 = nodeCoord(n1,:); p2 = nodeCoord(n2,:);
            h.allBars(e) = plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], '-', ...
                'Color', allC, 'LineWidth', allLW);
        end
    end

    % Highlight the selected bar
    n1 = barIJ(barID,1); n2 = barIJ(barID,2);
    p1 = nodeCoord(n1,:); p2 = nodeCoord(n2,:);
    h.hiBar = plot3(ax, [p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], '-', ...
        'Color', hiC, 'LineWidth', hiLW);

    % Label at midpoint
    if doLabel
        pm = 0.5*(p1+p2);
        h.txt = text(ax, pm(1), pm(2), pm(3), sprintf("  bar #%d", barID), ...
            'Color', hiC, 'FontWeight', 'bold');
    end

    if ~isempty(vw)
        view(ax, vw(1), vw(2));
    end

    hold(ax, 'off');
end