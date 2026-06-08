% Used to identify the critical points of the potential function and their classification

% --- Symbolic setup ---
syms x y real

% Gradient (first derivatives)
Fx = diff(F, x);
Fy = diff(F, y);

% Solve for critical points: grad F = 0
S = solve([Fx == 0, Fy == 0], [x, y], 'Real', true);

% Pack solutions into numeric array
crit_pts = double([S.x, S.y]);

% Hessian (second derivatives)
H = hessian(F, [x, y]);

% Helper for classification based on Hessian at a point
classifyPoint = @(Hp) classify_from_hessian(Hp);

% Evaluate and classify each critical point
fprintf('Critical points and classification for F(x,y):\n\n');
for k = 1:size(crit_pts, 1)
    xk = crit_pts(k,1);
    yk = crit_pts(k,2);

    % Evaluate Hessian and function value at the point
    Hk = double(subs(H, {x, y}, {xk, yk}));
    Fk = double(subs(F, {x, y}, {xk, yk}));

    % Classification using determinant and H(1,1)
    [label, D, a11] = classifyPoint(Hk);

    fprintf('Point %d: (x, y) = (%.6g, %.6g)\n', k, xk, yk);
    fprintf('  F(x,y) = %.6g\n', Fk);
    fprintf('  det(H) = %.6g,  H11 = %.6g\n', D, a11);
    fprintf('  Type:   %s\n\n', label);
end

% -------------------------
% Local function for classification
% -------------------------
function [label, D, a11] = classify_from_hessian(H)
    D   = det(H);
    a11 = H(1,1);

    if D > 0 && a11 > 0
        label = 'Local minimum';
    elseif D > 0 && a11 < 0
        label = 'Local maximum';
    elseif D < 0
        label = 'Saddle point';
    else
        label = 'Inconclusive (determinant ~ 0)';
    end
end