%% make date
clc
clear

fprintf('---------------------------PNP test---------------------------\n\n')

A = rand(3, 3);
t = rand(3, 1);
[R, ~] = qr(A);
P = [R, t];

npoints = 6;
X = rand(3, npoints);
X(4, :) = ones(1, npoints);

x = P * X;
for i = 1:npoints
    x(:, i) = x(:, i) ./ x(3, i);
end

errors = zeros(npoints, 4);

fprintf('number of points = %d\n', npoints)

%% add noise
std = 0.005;
noise = normrnd(0, std, 3, npoints);
X_noise = X;
X_noise(1:3, :) = X_noise(1:3, :) + noise;

P_norm = P ./ P(3, 4);

fprintf('noise std = %f\n\n', std)

%% PNP DTL
fprintf('---------------------------PNP DLT---------------------------\n\n')

P_est = solve_pnp_dlt(X_noise, x);
[K_est, R_est, t_est] = KRt_from_P(P_est);

angle_axis = rodrigues(R);
angle_axis_est = rodrigues(R_est);

fprintf('angle_axis = \n')
disp(angle_axis')
fprintf('angle_axis_est = \n')
disp(angle_axis_est')
fprintf('t = \n')
disp(t')
fprintf('t_est = \n')
disp(t_est')

x_est = [R_est, t_est] * X;
for i = 1:npoints
    x_est(:, i) = x_est(:, i) ./ x_est(3, i);
end
errors(:, 1) = sqrt(sum((x_est  - x).^2, 1))';
fprintf('reproject error = %f\n\n', mean(errors(:, 1)))

%% PNP P3P
fprintf('---------------------------PNP P3P---------------------------\n\n')

P_est = solve_pnp_p3p(X_noise(:, 1:3), x(:, 1:3));

if size(P_est, 1) ~= 0
    min_err = 1e10;
    min_i = 0;
    for i = 1:size(P_est, 2)
        [~, R_est, t_est] = KRt_from_P(P_est{i});

        angle_axis = rodrigues(R);
        angle_axis_est = rodrigues(R_est);
        err = sum(sum((angle_axis - angle_axis_est).^2)) + sum((t - t_est).^2);
        if err < min_err
            min_i = i;
            min_err = err;
        end
    end
    
    [K_est, R_est, t_est] = KRt_from_P(P_est{min_i});
    
    angle_axis = rodrigues(R);
    angle_axis_est = rodrigues(R_est);
    
    fprintf('angle_axis = \n')
    disp(angle_axis')
    fprintf('angle_axis_est = \n')
    disp(angle_axis_est')
    fprintf('t = \n')
    disp(t')
    fprintf('t_est = \n')
    disp(t_est')

    x_est = [R_est, t_est] * X;
    for i = 1:npoints
        x_est(:, i) = x_est(:, i) ./ x_est(3, i);
    end
    errors(:, 2) = sqrt(sum((x_est  - x).^2, 1))';
    fprintf('reproject error = %f\n\n', mean(errors(:, 2)))
else
    fprintf('p3p failed\n\n')
end

%% PNP EPnP
fprintf('---------------------------PNP EPnP---------------------------\n\n')

P_est = solve_pnp_epnp(X_noise, x);
[K_est, R_est, t_est] = KRt_from_P(P_est);

angle_axis = rodrigues(R);
angle_axis_est = rodrigues(R_est);

fprintf('angle_axis = \n')
disp(angle_axis')
fprintf('angle_axis_est = \n')
disp(angle_axis_est')
fprintf('t = \n')
disp(t')
fprintf('t_est = \n')
disp(t_est')

x_est = [R_est, t_est] * X;
for i = 1:npoints
    x_est(:, i) = x_est(:, i) ./ x_est(3, i);
end
errors(:, 3) = sqrt(sum((x_est  - x).^2, 1))';
fprintf('reproject error = %f\n\n', mean(errors(:, 3)))

%% PNP AP3P
fprintf('---------------------------PNP AP3P---------------------------\n\n')

P_est = solve_pnp_ap3p(X_noise(:, 1:3), x(:, 1:3));

if size(P_est, 1) ~= 0
    min_err = 1e10;
    min_i = 0;
    for i = 1:size(P_est, 2)
        [~, R_est, t_est] = KRt_from_P(P_est{i});

        angle_axis = rodrigues(R);
        angle_axis_est = rodrigues(R_est);
        err = sum(sum((angle_axis - angle_axis_est).^2)) + sum((t - t_est).^2);
        if err < min_err
            min_i = i;
            min_err = err;
        end
    end
    
    [K_est, R_est, t_est] = KRt_from_P(P_est{min_i});
    
    angle_axis = rodrigues(R);
    angle_axis_est = rodrigues(R_est);
    
    fprintf('angle_axis = \n')
    disp(angle_axis')
    fprintf('angle_axis_est = \n')
    disp(angle_axis_est')
    fprintf('t = \n')
    disp(t')
    fprintf('t_est = \n')
    disp(t_est')

    x_est = [R_est, t_est] * X;
    for i = 1:npoints
        x_est(:, i) = x_est(:, i) ./ x_est(3, i);
    end
    errors(:, 4) = sqrt(sum((x_est  - x).^2, 1))';
    fprintf('reproject error = %f\n\n', mean(errors(:, 4)))
else
    fprintf('ap3p failed\n\n')
end

boxplot(errors)
