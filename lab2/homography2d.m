function H = homography2d(x1, x2)
    % H = homography2d(x1, x2)
    %
    % x1: set of points in plane 1. size(x1) = [ncoord, npoints]
    % x2: set of points in plane 2. size(x2) = [ncoord, npoints]
    %
    % This function takes two sets of corresponding points between two
    % projective planes P^2, and compute the homography tha relates them.
    % If the system does not have an exact solution, the least squares
    % solution will be found.
    
    % Number of points:
    npoints = size(x1,2);
    
    % Checkings for consistency:
    if(npoints < 4)
        error('Not enough points for determining the homography.')
    end
    if(npoints ~= size(x2,2))
        error('Different number of points in each plane.')
    end
    if(size(x1,1) ~= 3 || size(x2,1) ~= 3)
        error('Points must have three coordinates.')
    end
    
    % Initialize A:
    A = zeros(2 * npoints, 9);
    
    % Loop over each pair of correspondences:
    for i = 1:npoints
        A((2*i-1):(2*i), :) = [0, 0, 0           , -x2(3,i) * x1(:,i)', x2(2,i) * x1(:,i)';
                               x2(3,i) * x1(:,i)', 0, 0, 0            , -x2(1,i) * x1(:,i)'];
    end
    
    % Singula value decomposition:
    [~, ~, V] = svd(A);
    
    % Least squares solutions (exact if the system is not over-determined):
    h = V(:,9);

    % Build the matrix of the proyective transformation:
    H = [h(1), h(2), h(3);
         h(4), h(5), h(6);
         h(7), h(8), h(9)];

    return
end