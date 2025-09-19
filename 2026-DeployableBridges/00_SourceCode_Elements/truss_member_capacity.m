function result = truss_member_capacity(N, Ag, An, Fy, Fu, E, K, L, r, phi)
% 计算三维桁架每根杆的承载力（LRFD）
% N  : [nBars x 1] 轴力（+拉，-压），与 internal_force 一致的单位
% Ag : 毛截面面积
% An : 净截面面积（若无孔洞则设为 Ag）
% Fy : 屈服强度
% Fu : 抗拉强度
% E  : 弹性模量
% K  : 计算长度系数（无特别约束可取 1）
% L  : 杆长
% r  : 回转半径（取弱轴，若实心圆可用 d/4）
% phi: 结构体，含字段 .t_y, .t_u, .c（默认 0.90/0.75/0.90）
%
% 返回：
% result.R        : 设计承载力（与 N 同单位）
% result.util     : 利用系数 |N|/R
% result.mode     : 控制破坏模式（cellstr）
% result.slender  : KL/r
% result.flags.*  : 细长比经验限值通过与否
% result.details.*: 中间量（便于调试/查验）

    % ---- 默认的 phi ----
    if nargin < 10 || isempty(phi)
        phi = struct('t_y',0.90, 't_u',0.75, 'c',0.90);
    end

    % ---- 转列向量并做标量扩展 ----
    makecol = @(x) double(x(:));
    N  = makecol(N);   n = numel(N);

    Ag = makecol(default_to(Ag,  1));
    An = makecol(default_to(An, Ag));
    Fy = makecol(default_to(Fy,  1));
    Fu = makecol(default_to(Fu, Fy));
    E  = makecol(default_to(E,   2e11));
    K  = makecol(default_to(K,   1));
    L  = makecol(default_to(L,   1));
    r  = makecol(default_to(r,   1));

    % 扩展标量到 nBars
    Ag = expand_to(Ag,n); An = expand_to(An,n);
    Fy = expand_to(Fy,n); Fu = expand_to(Fu,n);
    E  = expand_to(E,n);  K  = expand_to(K,n);
    L  = expand_to(L,n);  r  = expand_to(r,n);

    % 数值健壮性（避免 0 长度/半径导致的除零）
    K(K<=0) = 1;
    L(L<=0) = eps;
    r(r<=0) = eps;

    % ---- 受拉承载力：毛截面屈服 vs 净截面断裂 ----
    Pn_y      = Fy .* Ag;                 % 毛截面屈服
    Pn_u      = Fu .* An;                 % 净截面断裂
    phiPn_t_y = phi.t_y .* Pn_y;
    phiPn_t_u = phi.t_u .* Pn_u;
    phiPn_t   = min(phiPn_t_y, phiPn_t_u);  % 逐元素最小值

    % 哪个受拉极限状态控制
    t_yield_ctrl = (phiPn_t_y <= phiPn_t_u);

    % ---- 受压承载力：AISC 压杆分段公式 ----
    slender  = (K .* L) ./ r;                          % KL/r
    Fe       = (pi.^2 .* E) ./ max(slender.^2, eps);   % 欧拉应力
    lambda_c = sqrt(Fy ./ max(Fe, eps));               % 无量纲细长度

    Fcr = zeros(n,1);
    inelastic = (lambda_c <= 1.5);
    Fcr(inelastic)  = (0.658 .^ (lambda_c(inelastic).^2)) .* Fy(inelastic);
    Fcr(~inelastic) = 0.877 .* Fe(~inelastic);

    phiPn_c = phi.c .* Fcr .* Ag;

    % ---- 按内力正负选取承载力 ----
    isT = (N >= 0); isC = ~isT;
    R = zeros(n,1);
    R(isT) = phiPn_t(isT);
    R(isC) = phiPn_c(isC);

    % ---- 模式标记 ----
    mode = repmat({''}, n, 1);
    mode(isT &  t_yield_ctrl) = {'tension_yield'};
    mode(isT & ~t_yield_ctrl) = {'tension_fracture'};
    mode(isC)                 = {'compression_buckling'};

    % ---- 利用系数 ----
    R_safe = max(R, realmin);     % 防止除零
    util   = abs(N) ./ R_safe;

    % ---- 经验细长比限值（服务性，不是强度）----
    slender_ok_comp = true(n,1);
    slender_ok_tens = true(n,1);
    slender_ok_comp(isC) = slender(isC) <= 200;
    slender_ok_tens(isT) = slender(isT) <= 300;

    % ---- 输出 ----
    result.R        = R;
    result.util     = util;
    result.mode     = mode;
    result.slender  = slender;
    result.flags.slender_ok_comp = slender_ok_comp;
    result.flags.slender_ok_tens = slender_ok_tens;
    result.details.phiPn_t_y = phiPn_t_y;
    result.details.phiPn_t_u = phiPn_t_u;
    result.details.phiPn_c   = phiPn_c;
    result.details.Fcr       = Fcr;
end

% ------- 小工具函数（同文件内的子函数） -------
function x = default_to(x, default_val)
    if nargin==0 || isempty(x), x = default_val; end
end

function y = expand_to(x, n)
    if numel(x)==1
        y = repmat(x, n, 1);
    elseif numel(x)==n
        y = x;
    else
        error('Input size mismatch: expected scalar or length %d, got %d.', n, numel(x));
    end
end
