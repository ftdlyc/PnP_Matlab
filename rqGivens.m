function [ R, Q ] = rqGivens( A )
% rqGivens Calculates RQ decomposition of A = RQ (3x3)
%
% Syntax:
%   [R, Q] = rqGivens(A);
%
% Input:
%   A - 3-by-3 matrix of rank 3
%
% Output:
%   R - Upper triangular matrix (3-by-3)
%   Q - Orthogonal matrix (3-by-3)
%
% Description:
%   This function calculates the 3-dimensional RQ decomposition of A using 
%   Givens rotations (equal to Euler rotations) Gx, Gy Gz:
%
%   Gx = [  1   0   0;
%           0   c   -s;
%           0   s   c];
%     
%   Gy = [  c   0   s;
%           0   1   0;
%           -s  0   c];
%
%   Gz = [  c   -s  0;
%           s   c   0;
%           0   0   1];
%
%   Ax = A * Gx to set Ax(3,2) to zero.
%   Axy = Ax * Gy to set Axy(3,1) to zero.
%   R = Axyz = Axy * Gz to set Axyz(2,1) to zero.
%
%   R = A * Gx * Gy * Gz 
%   -> R * Gz' * Gy' * Gx' = A
%   -> Q = Gz' * Gy' * Gx'
%
%   See also: 
%   - https://en.wikipedia.org/wiki/Givens_rotation#Dimension_3
%   - Hartley, Zisserman - Multiple View Geometry in Computer Vision
%     http://www.amazon.com/dp/0521540518 (Appendix 4, A4.1.1, page 579)     
%
% Author: Lars Meinel
% Email: lars.meinel@etit.tu-chemnitz.de
% Date: July 2015; Last revision: 2015-07-10
%% 1st - Set element 32 to zero
if A(3,2) == 0
    Ax = A;
    Gx = eye(3);
else
    r32 = sqrt(A(3,3)*A(3,3) + A(3,2)*A(3,2));
    c32 = A(3,3) / r32;
    s32 = -A(3,2) / r32;
    G32 = [ 1    0    0;
            0    c32  -s32;
            0    s32  c32   ];
    Gx = G32;   
    Ax = A * Gx;
end
%% 2nd - Set element 31 to zero
if A(3,1) == 0
    Axy = Ax;
    Gy = eye(3);
else
    r31 = sqrt(Ax(3,3)*Ax(3,3) + Ax(3,1)*Ax(3,1));
    c31 = Ax(3,3) / r31;
    s31 = Ax(3,1) / r31;
    G31 = [ c31     0   s31;
            0       1   0;
            -s31    0   c31];
    Gy = G31;   
    Axy = Ax * Gy;
end
%% 3rd - Set element 21 to zero
if A(2,1) == 0
    Axyz = Axy;
    Gz = eye(3);
else
    r21 = sqrt(Axy(2,2)*Axy(2,2) + Axy(2,1)*Axy(2,1));
    c21 = Axy(2,2) / r21;
    s21 = -Axy(2,1) / r21;
    G21 =    [  c21     -s21    0;
                s21     c21     0;
                0       0       1];
    Gz = G21;
    Axyz = Axy * Gz;
end
R = Axyz;
Q = Gz' * Gy' * Gx';
end