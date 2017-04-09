function [Pproj, Xproj] = factorization_method(x)
% computes a projective reconstruction with the factorization method of 
% Sturm and Triggs '1996
%
% Arguments:
%     - x: cell array with as many elements as the number of cameras. Each
%     element is a 3 x Npoints matrix, with the homogeneous coordinates of
%     the 2D points in the correspondent camera.
%
% Returns:
%     - Pproj: 3*Ncameras x 4 matrix with the projective matrices of all
%     cameras (each 3 rows is one camera matrix).
%     - Xproj: 4*Npoints x 1 array with the homogeneous coordinates of the
%     3D points.

    % Some parameters:
    tol_reescale = 0.1;
    tol_d = 0.1;

    % Some numbers:
    Ncam = length(x);
    Npoints = size(x{1},2);

    % Normalize points:
    % TODO

    % Initialize lambdas:
    lambda = ones(Ncam, Npoints);
    
    % Initialize mean distance between original and projected 2D poitns:
    d = Inf;
    
    % Loop.
    iterate = 1;
    while(iterate)
        % Reescale matrix lambda:
        reescale = 1;
        while(reescale)
            lambda_old = lambda;
            % Reescale rows:
            for i = 1:Ncam
                lambda(i,:) = lambda(i,:) / norm(lambda(i,:));
            end
            % Reescale columns:
            for j = 1:Npoints
                lambda(:,j) = lambda(:,j) / norm(lambda(:,j));
            end
            % Stop condition:
            if(norm(lambda - lambda_old) < tol_reescale)
                reescale = 0;
            end
        end
            
        % Build measurement matrix:
        M = zeros(3*Ncams, Npoints);
        for i = 1:Ncams
            M(3*i-2,:) = lambda(i,:) .* x{i}(1,:);
            M(3*i-1,:) = lambda(i,:) .* x{i}(2,:);
            M(3*i,:)   = lambda(i,:) .* x{i}(3,:);
        end
        
        % Singular Value Decomposition of M:
        [U,D,V] = svd(M);
        
        % Camera projection matrices and homogenous 3D points:
        Pproj = U * D(:,1:4);
        Xproj = V(:,1:4)';
        
        % Compute distance between original 2D points and projected ones:
        d_old = d;
        d = 0;
        for i = 1:Ncams
            x_proj_2d = Pproj((3*i-2):(3*i),:) * Xproj;
            x_proj_euc = euclid(x_proj_2d);
            for j = 1:Npoints
                d = d + norm(x{i}(:,j) - x_proj_euc(:,j))^2;
            end
        end
        % Divide to get the mean distance (although it does not affect the
        % result at all):
        d = d / (Npoints * Ncams);
        
        % Exit condition:
        if(abs(d_old - d) / d < tol_d)
            iterate = 0;
        end
    end
    
    % Unnormalize:
    % TODO
    
    % Triangulate and resection ???

    return
    
end