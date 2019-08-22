% Perspective-n-Point
% AP3P Method 
% Reference: An Efficient Algebraic Solution to the Perspective-Three-Point Problem
%
% by ftdlyc
%
% Input
% X: [3 x n] or [4 x n] 3D points
% x: [2 x n] or [3 x n] 2D points
%
% Output
% P: n x [3 x 4] cell of Camera Projection Martix
%
function P = solve_pnp_p3p(X, x)
P = cell(0, 0);

[row_X, col_X] = size(X);
[row_x, col_x] = size(x);
if row_X ~= 3 && row_X ~= 4
    fprintf('mat X must be [3 x n] or [4 x n]\n')
    return
end
if row_x ~= 2 && row_x ~= 3
    fprintf('mat x must be [2 x n] or [3 x n]\n')
    return
end
if col_X ~= col_x
    fprintf('col(x) no equal to col(X)\n')
    return
end
n = col_X;
if n ~= 3
    fprintf('col(x)/col(X) must be 3\n')
    return
end

if row_X == 4
    X = X(1:3, :) ./ X(4, :);
end
if row_x == 2
    x(3, :) = ones(1, col_x);
end

p1 = X(:, 1);
p2 = X(:, 2);
p3 = X(:, 3);
b1 = x(:, 1);
b2 = x(:, 2);
b3 = x(:, 3);

%% compute k1, k3
k1 = (p1 - p2) / norm(p1 - p2);
k3 = cross(b1, b2);
nk3 = norm(k3);
k3 = k3 / nk3;

%% compute u1, u2, v1, v2
u1 = p1 - p3;
u2 = p2 - p3;
v1 = cross(b1, b3);
v2 = cross(b2, b3);

%% compute delta, nl
nl = cross(u1, k1);
delta = norm(nl);
nl = nl / delta;

%% compute f
k3b3 = k3' * b3;
b1b2 = b1' * b2;

f11 = delta * k3b3;
f21 = delta * b1b2 * k3b3;
f22 = delta * k3b3 * nk3;
f13 = delta * v1' * k3;
f23 = delta * v2' * k3;
f24 = u2' * k1 * k3b3 * nk3;
f15 = - u1' * k1 * k3b3;
f25 = - u2' * k1 * b1b2 * k3b3;

%% compute quartic polynomial coeffs
g1 = f13 * f22;
g2 = f13 * f25 - f15 * f23;
g3 = f11 * f23 - f13 * f21;
g4 = -f13 * f24;
g5 = f11 * f22;
g6 = f11 * f25 - f15 * f21;
g7 = -f15 * f24;

coeff4 = g5 * g5 + g1 * g1 + g3 * g3;
coeff3 = 2 * (g5 * g6 + g1 * g2 + g3 * g4);
coeff2 = g6 * g6 + 2 * g5 * g7 + g2 * g2 + g4 * g4 - g1 * g1 - g3 * g3;
coeff1 = 2 * (g6 * g7 - g1 * g2 - g3 * g4);
coeff0 = g7 * g7 - g2 * g2 - g4 * g4;

%% solve quartic polynomial and refine
s = roots([coeff4, coeff3, coeff2, coeff1, coeff0]);
real_roots = zeros(4, 1);

niters = 2;
for i = 1:4
    root = real(s(i));
    for j = 1:2
      error = coeff0 + root * (coeff1 + root * (coeff2 + root * (coeff3 + root * coeff4)));
      derivative = coeff1 + root * (2 * coeff2 + root * (3 * coeff3 + 4 * root * coeff4));
      root = root - error / derivative;
    end
    real_roots(i) = root;
end

%% compute C and pc
Ck1nl = [k1 nl cross(k1, nl)];
Cb1k3T = [b1 k3 cross(b1, k3)]';

n = 1;
for i = 1:4
    ctheta1 = real_roots(i);
    if abs(ctheta1) > 1
        continue
    end
    
    stheta1 = sign(k3b3) * sqrt(1 - ctheta1 * ctheta1);
    ntheta3 = stheta1 / (g5 * ctheta1 * ctheta1 + g6 * ctheta1 + g7);
    ctheta3 = ntheta3 * (g1 * ctheta1 + g2);
    stheta3 = ntheta3 * (g3 * ctheta1 + g4);
    
    C = Ck1nl ...
        * [1     0       0   ; ...
           0  ctheta1 stheta1; ...
           0 -stheta1 ctheta1] ...
        * [ctheta3 0 -stheta3; ...
              0    1     0   ; ...
           stheta3 0  ctheta3] ...
        * Cb1k3T;
    
    R = C';
    t = delta * stheta1 / k3b3 * b3 - R * p3;
    P{n} = [R, t];
    n = n + 1;
end
end
