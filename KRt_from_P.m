% Extract K, R, t from Projection Martix
%
% by ftdlyc
%
% Input
% P: [3 x 4] Camera Projection Martix
%
% Output
% K: [3 x 3] Camera Intrinsic Matrix
% R: [3 x 3] Rotation Matrix
% t: [1 x 3] Translation Vector
%
function [K, R, t] = KRt_from_P(P)
%% QR decomposition
[K, R] = rqGivens(P(1:3, 1:3));

%% ensure that the diagonal is positive
if K(3, 3) < 0
    K = -K;
    R = -R;
end
if K(2, 2) < 0
    S = [1  0  0 
         0 -1  0
         0  0  1];
    K = K * S;
    R = S * R;
end
if K(1, 1) < 0
    S = [-1  0  0 
          0  1  0
          0  0  1];
    K = K * S;
    R = S * R;
end

%% ensure R determinant == 1
t = linsolve(K, P(:, 4));

if det(R) < 0
    R = -R;
    t = -t;
end

K = K ./ K(3, 3);

end
