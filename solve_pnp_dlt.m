% Perspective-n-Point
% Direct Linear Transform
%
% by ftdlyc
%
% Input
% X: [3 x n] or [4 x n] 3D points
% x: [2 x n] or [3 x n] 2D points
%
% Output
% P: [3 x 4] Camera Projection Martix
%
function P = solve_pnp_dlt(X, x)
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
if n < 6
    fprintf('col(x)/col(X) must greather than 6\n')
    return
end

if row_X == 3
    X(4, :) = ones(1, col_X);
end

%% build action matrix
A = zeros(n * 2, 12);
for i = 1:n
    pt3d = X(:, i)';
    A(2 * i - 1, :) = [pt3d,         zeros(1, 4),  -x(1, i) .* pt3d];
    A(2 * i, :)     = [zeros(1, 4),  pt3d,         -x(2, i) .* pt3d];
end

%% Solve P
[~, D, V] = svd(A);
radio = D(11, 11) / D(1, 1);
if radio > 1e-5
    P = V(:, 12);
    P = reshape(P, 4, 3)';
    P = P ./ P(3, 4);
else
    P = [eye(3) zeros(3, 1)];
end

end
