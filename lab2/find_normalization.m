function [xnorm, Hnorm] = find_normalization(x)

    % First we translate to have the centroid in the origin, and the we
    % scale to ensure the average distance to the origin is sqrt(2).

    tx = mean(x(1,:) ./ x(3,:));
    ty = mean(x(2,:) ./ x(3,:));

    % Transformation for centering the data (translation):
    Ht = [1, 0, -tx;
         0, 1, -ty;
         0, 0, 1];

    % Translated data:
    xt = Ht * x;

    % Average distance to the origin:
    d_xt = avg_dist(xt);

    % Solve the equation for the scaling factor:
    s = sqrt(2) / d_xt;

    % Scaling transformation:
    Hs = [s, 0, 0;
         0, s, 0;
         0, 0, 1];

    % Combined transformation:
    Hnorm = Hs * Ht;
    
    % Normalized data:
    xnorm = Hnorm * x; % should be the same as Hs * st

    return
end