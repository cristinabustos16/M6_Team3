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

% lpmayos note: theory lecture 8 slide 50

    % Some parameters:
    tol_reescale = 0.1;
    tol_d = 0.1;

    % Some numbers:
    Ncams = length(x);
    Npoints = size(x{1},2);

    % 2. Normalize the set of points in each image (similarity transf. Hs).
    x_hat = cell(1,Ncams);
    Harray = cell(1,Ncams);
    for i = 1:Ncams
        [x_hat{i}, Harray{i}] = normalise2dpts(x{i});
    end

    % 3. Initialize all lambda_ij (= 1 or better initialization).

    lambda = ones(Ncams, Npoints);

    % Sturm initialization (NOTE: initialiing to 1 does not work!)
    
    % camera 1
    F1 = fundamental_matrix(x{1}, x{1});
    [U, D, V] = svd(F1);
    e1 = V(:,3) / V(3,3);
        
    % camera 2
    F2 = fundamental_matrix(x{2}, x{1});
    [U, D, V] = svd(F2);
    e2 = V(:,3) / V(3,3);
        
    for j=1:size(x{1},2)
        num = x{1}(:, j)' * F1 * cross(e1, x{1}(:,j));
        denom = norm(cross(e1, x{1}(:,j))) .^ 2 * lambda(1, j);
        lambda(1,j) = num / denom;
    end
    for j=1:size(x{2},2)
        num = x{1}(:, j)' * F2 * cross(e2, x{2}(:,j));
        denom = norm(cross(e2, x{2}(:,j))) .^ 2 * lambda(1, j);
        lambda(2,j) = num / denom;
    end

    % Initialize mean distance between original and projected 2D poitns:
    d = Inf;
    
    % Loop.
    iterate = 1;
    while(iterate)

        % 4. Alternate rescaling the rows of the depth matrix (formed by lambda_ij)
        % to have unit norm and the columns of depth matrix to have unit norm until
        % depth matrix stops changing significantly (usually two loops).
        % Reescale matrix lambda:
        reescale = 1;
        while(reescale)
            lambda_old = lambda;
            % Reescale rows:
            for i = 1:Ncams
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
            
        % 5. Build the measurement matrix M
        M = zeros(3*Ncams, Npoints);
        for i = 1:Ncams
            M(3*i-2,:) = lambda(i,:) .* x_hat{i}(1,:);
            M(3*i-1,:) = lambda(i,:) .* x_hat{i}(2,:);
            M(3*i,:)   = lambda(i,:) .* x_hat{i}(3,:);
        end
        
        % 6. Determine the SVD of M = UDVT
        [U,D,V] = svd(M);
        
        % 7. Let PM = UD4 and XM =V4T
        % Camera projection matrices and homogenous 3D points:
        Pproj_hat = U * D(:,1:4);
        Xproj = V(:,1:4)';
        
        % 8. If sum_i sum_j d(xji , P^i X_j)2 converges then stop; 
        % Otherwise let lambda_ij = (P^iX_j)_3 and go to Step 4.
        % Compute distance between original 2D points and projected ones:
        d_old = d;
        d = 0;
        for i = 1:Ncams
            x_proj_2d = Pproj_hat((3*i-2):(3*i),:) * Xproj;
            x_proj_euc = euclid(x_proj_2d);
            x_orig_euc = euclid(x_hat{i});
            for j = 1:Npoints
                d = d + norm(x_orig_euc(:,j) - x_proj_euc(:,j))^2;
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
    
    % 9. Unnormalize the camera matrices
    for i = 1:Ncams
        Pproj((3*i-2):(3*i),:) = inv(Harray{i}) * Pproj_hat((3*i-2):(3*i),:);
    end
    
    % Triangulate and resection ???

    return
    
end