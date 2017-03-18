function d = avg_dist(x)

    % Compute the average distance to the origin of a set of points.
    d = mean(sqrt((x(1,:).^2 + x(2,:).^2) / x(3,:).^2));
    
    return
end