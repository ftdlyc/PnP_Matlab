% Perspective-n-Point
% EPnP Method 
% Reference: EPnP An Accurate O(n) Solution to the PnP Problem
%
% by ftdlyc
%
% Input
% X: [3 x n] or [4 x n] 3D points
% x: [2 x n] or [3 x n] 2D points
%
% Output
% P: [3 x 4] cell of Camera Projection Martix
%
function P = solve_pnp_epnp(X, x)
err = 1e10;
R = eye(3);
t = zeros(3, 1);
P = [R t];

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
if n <= 3
    fprintf('col(x)/col(X) must greather than 3\n')
    return
end

if row_X == 4
    X = X(1:3, :) ./ X(4, :);
end
if row_x == 3
    x = x(1:2, :) ./ x(3, :);
end

%% choose control point
C = zeros(3, 4);
C(:, 1) = sum(X, 2) ./ n;
A = X - C(:, 1);
[v, lambda] = eig(A * A');
for i = 1:3
    C(:, i + 1) = C(:, 1) + sqrt(lambda(i, i) / n) * v(:, i);
end

%% compute barycentric coordinates
% 4 x n matrix
% [X; 1] = [C; 1] * alpha
alpha = pinv([C; ones(1, 4)]) * [X; ones(1, n)];

%% compute V
fx = 1;
fy = 1;
cx = 0;
cy = 0;
M = zeros(2 * n, 12);
for i = 1:n
    M(2 * i - 1, :) = [alpha(1, i) * fx, 0, alpha(1, i) * (cx - x(1, i)) ...
                       alpha(2, i) * fx, 0, alpha(2, i) * (cx - x(1, i)) ...
                       alpha(3, i) * fx, 0, alpha(3, i) * (cx - x(1, i)) ...
                       alpha(4, i) * fx, 0, alpha(4, i) * (cx - x(1, i))];
    M(2 * i, :) = [0, alpha(1, i) * fy, alpha(1, i) * (cy - x(2, i)) ...
                   0, alpha(2, i) * fy, alpha(2, i) * (cy - x(2, i)) ...
                   0, alpha(3, i) * fy, alpha(3, i) * (cy - x(2, i)) ...
                   0, alpha(4, i) * fy, alpha(4, i) * (cy - x(2, i))];
end
[~, ~, V] = svd(M);

%% compute L_6x10 and rho
% betas10 = [b11 b12 b22 b13 b23 b33 b14 b24 b34 b44]
% L = [l12'; l13'; l14'; l23'; l24'; l34';]
v = zeros(12, 4);
v(:, 1) = V(:, 12);
v(:, 2) = V(:, 11);
v(:, 3) = V(:, 10);
v(:, 4) = V(: , 9);

s = cell(4, 6);
for i = 1:4
    vtmp1 = v(1:3, i);
    vtmp2 = v(4:6, i);
    vtmp3 = v(7:9, i);
    vtmp4 = v(10:12, i);
    s{i, 1} = vtmp1 - vtmp2;
    s{i, 2} = vtmp1 - vtmp3;
    s{i, 3} = vtmp1 - vtmp4;
    s{i, 4} = vtmp2 - vtmp3;
    s{i, 5} = vtmp2 - vtmp4;
    s{i, 6} = vtmp3 - vtmp4;
end

L = zeros(6, 10);
for i = 1:6
    L(i, :) = [    s{1, i}' * s{1, i} ...
               2 * s{1, i}' * s{2, i} ...
                   s{2, i}' * s{2, i} ...
               2 * s{1, i}' * s{3, i} ...
               2 * s{2, i}' * s{3, i} ...
                   s{3, i}' * s{3, i} ...
               2 * s{1, i}' * s{4, i} ...
               2 * s{2, i}' * s{4, i} ...
               2 * s{3, i}' * s{4, i} ...
                   s{4, i}' * s{4, i}];
end

rho = zeros(6 , 1);
rho(1) = sum((C(:, 1) - C(:, 2)).^2);
rho(2) = sum((C(:, 1) - C(:, 3)).^2);
rho(3) = sum((C(:, 1) - C(:, 4)).^2);
rho(4) = sum((C(:, 2) - C(:, 3)).^2);
rho(5) = sum((C(:, 2) - C(:, 4)).^2);
rho(6) = sum((C(:, 3) - C(:, 4)).^2);

%% compute beta (N = 1)
% find betas
betas = zeros(4, 1);
b1 = linsolve(L(:, 1), rho);
betas(1) = sqrt(abs(b1));

betas = gauss_newton(L, rho, betas);

Xc = compute_Xc(alpha, betas, v);
[R1, t1] = compute_Rt(X, Xc);
err1 = reproject_error(X, x, R1, t1);

if err1 < err
    err = err1;
    R = R1;
    t = t1;
end

%% compute beta (N = 2)
% betas10 = [b11 b12 b22 b13 b23 b33 b14 b24 b34 b44]
% b3      = [b11 b12 b22]
% find beta
betas = zeros(4, 1);
b3 = linsolve(L(:, 1:3), rho);
betas(1) = sqrt(abs(b3(1)));
betas(2) = sqrt(abs(b3(3))) * sign(b3(2)) * sign(b3(1));

betas = gauss_newton(L, rho, betas);

Xc = compute_Xc(alpha, betas, v);
[R2, t2] = compute_Rt(X, Xc);

err2 = reproject_error(X, x, R2, t2);
if err2 < err
    err = err2;
    R = R2;
    t = t2;
end

%% compute beta (N = 3)
% betas10 = [b11 b12 b22 b13 b23 b33 b14 b24 b34 b44]
% b6      = [b11 b12 b22 b13 b23 b33]
% find beta
betas = zeros(4, 1);
b6 = linsolve(L(:, 1:6), rho);
betas(1) = sqrt(abs(b6(1)));
betas(2) = sqrt(abs(b6(3))) * sign(b6(2)) * sign(b6(1));
betas(3) = sqrt(abs(b6(6))) * sign(b6(4)) * sign(b6(1));

betas = gauss_newton(L, rho, betas);

Xc = compute_Xc(alpha, betas, v);
[R3, t3] = compute_Rt(X, Xc);

err3 = reproject_error(X, x, R3, t3);
if err3 < err
    err = err3;
    R = R3;
    t = t3;
end

%% compute beta (N = 4)
% betas10 = [b11 b12 b22 b13 b23 b33 b14 b24 b34 b44]
% b4      = [b11 b12     b13         b14]
% find beta
betas = zeros(4, 1);
b4 = linsolve([L(:, 1), L(:, 2), L(:, 4), L(:, 7)], rho);
betas(1) = sqrt(abs(b4(1)));
betas(2) = b4(2) / betas(1);
betas(3) = b4(3) / betas(1);
betas(4) = b4(4) / betas(1);

betas = gauss_newton(L, rho, betas);

Xc = compute_Xc(alpha, betas, v);
[R4, t4] = compute_Rt(X, Xc);

err4 = reproject_error(X, x, R4, t4);
if err4 < err
    err = err4;
    R = R4;
    t = t4;
end

%% return
P = [R t];
P = P ./ P(3, 4);
end

function betas = gauss_newton(L, rho, betas)
n_iter = 5;
for i = 1:n_iter
    % 6 x 10
    % err 12, 13, 14, 23, 24, 34
    % betas 11 12 22 13 23 33 14 24 34 44
    J = zeros(6, 4);
    for j = 1:6
        J(j, 1) = 2 * L(j, 1) * betas(1) +     L(j, 2) * betas(2) +     L(j, 4) * betas(3) +     L(j, 7) * betas(4);
        J(j, 2) =     L(j, 2) * betas(1) + 2 * L(j, 3) * betas(2) +     L(j, 5) * betas(3) +     L(j, 8) * betas(4);
        J(j, 3) =     L(j, 4) * betas(1) +     L(j, 5) * betas(2) + 2 * L(j, 6) * betas(3) +     L(j, 9) * betas(4);
        J(j, 4) =     L(j, 7) * betas(1) +     L(j, 8) * betas(2) +     L(j, 9) * betas(3) + 2 * L(j, 10) * betas(4);
    end
    
    % 6 x 1
    % err 12, 13, 14, 23, 24, 34
    r = zeros(6, 1);
    for j = 1:6
        r(j) = rho(j) ...
             - L(j, 1) * betas(1) * betas(1) ...
             - L(j, 2) * betas(1) * betas(2) ...
             - L(j, 3) * betas(2) * betas(2) ...
             - L(j, 4) * betas(1) * betas(3) ...
             - L(j, 5) * betas(2) * betas(3) ...
             - L(j, 6) * betas(3) * betas(3) ...
             - L(j, 7) * betas(1) * betas(4) ...
             - L(j, 8) * betas(2) * betas(4) ...
             - L(j, 9) * betas(3) * betas(4) ...
             - L(j, 10) * betas(4) * betas(4);
    end 
    
    A = J' * J;
    b = J' * r;
    dbetas = linsolve(A, b);
    
    betas = betas + dbetas;
end
end

function Xc = compute_Xc(alpha, betas, v)
x = betas(1) * v(:, 1) + betas(2) * v(:, 2) + betas(3) * v(:, 3) + betas(4) * v(:, 4);
C = [x(1:3, 1), x(4:6, 1), x(7:9, 1), x(10:12, 1)];
Xc = C * alpha;
end

function [R, t] = compute_Rt(Xw, Xc)
mean_Xw = mean(Xw, 2);
mean_Xc = mean(Xc, 2);

new_Xw = Xw - mean_Xw;
new_Xc = Xc - mean_Xc;

W = zeros(3, 3);
for j = 1:size(Xw, 2)
    W = W + new_Xc(:, j) * new_Xw(:, j)'; 
end
    
[U, ~, V] = svd(W);
R = U * V';
t = mean_Xc - R * mean_Xw;
end

function err = reproject_error(X, x, R, t)
x2 = R * X + t;
for i = 1:size(X, 2)
    x2(:, i) = x2(:, i) ./ x2(3, i);
end
x2 = x2(1:2, :);
err = mean(sqrt(sum((x2 - x).^2, 1)));
end
