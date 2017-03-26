function [F, inliers] = ransac_fundamental_matrix(p1, p2, th)

% ransac
[Ncoords, Npoints] = size(p1);
it = 0;
best_inliers = [];
% probability that at least one random sample set is free of outliers
p = 0.999; 
max_it = 1000;
while it < max_it
    
    points = randomsample(Npoints, 8);
    F = fundamental_matrix(p1(:,points), p2(:,points));

    inliers = compute_inliers(F, p1, p2, th);
    
    % test if it is the best model so far
    if length(inliers) > length(best_inliers)
        best_inliers = inliers;
    end    
    
    % update estimate of max_it (the number of trials) to ensure we pick, 
    % with probability p, an initial data set with no outliers
    fracinliers =  length(inliers)/Npoints;
    pNoOutliers = 1 -  fracinliers^4;
    pNoOutliers = max(eps, pNoOutliers);  % avoid division by -Inf
    pNoOutliers = min(1-eps, pNoOutliers);% avoid division by 0
    max_it = log(1-p)/log(pNoOutliers)
    
    it = it + 1
end

% compute F from the best inliers
F = fundamental_matrix(p1(:,best_inliers), p2(:,best_inliers));
inliers = best_inliers;


% The inliers are obtained with a threshold on the first order 
% approximation of the geometric error: Sampson distance
function idx_inliers = compute_inliers(F, p1, p2, th)

     % transformed points (in both directions)
    Fx = F * p1;  % naming convention to follow slide 3 of lab3.pdf
    Ftxp = F' * p2;  % naming convention to follow slide 3 of lab3.pdf
    
    % http://cariparo.dei.unipd.it/documents/corso_psc_07-08/sensorfusion/matlab/filemcv/sampson.m/view
    p2tFx_size = size(p2,2);

    p2tFx = zeros(1,p2tFx_size);
    for i = 1:p2tFx_size
        p2tFx(i) = p2(:,i)'*Fx(:,i);
    end

    
    % compute the Sampson distance
    d = (p2tFx.^2) ./ (Fx(1,:).^2 + Fx(2,:).^2 + Ftxp(1,:).^2 + Ftxp(2,:).^2);
    idx_inliers = find(d < th.^2);



function item = randomsample(npts, n)
    a = [1:npts]; 
    item = zeros(1,n);    
    for i = 1:n
      % Generate random value in the appropriate range 
      r = ceil((npts-i+1).*rand);
      item(i) = a(r);       % Select the rth element from the list
      a(r)    = a(end-i+1); % Overwrite selected element
    end                       % ... and repeat