function error = gs_errfunction( P0, Xobs )

% get the homography
H = reshape(P0(1:9), [3,3]);

points_lenght = size(Xobs,1) / 2;

x = Xobs(1:points_lenght);
x = reshape(x, [2,size(x,1)/2]);

xp = Xobs(points_lenght+1:end);
xp = reshape(xp, [2,size(xp,1)/2]);

%Get x_hat
x_hat = P0(10:end); % the first 9 elements correspond to the homography
lenght_x_hat = size(x_hat,1)/2;
x_hat = reshape(x_hat, [2,lenght_x_hat]);
x_hat = [x_hat; ones(1,lenght_x_hat)];


% Calculate Homography * x_hat and the error between estimates and observations
x_hat_p = H * x_hat;

error1 = abs(x - euclid(x_hat));
error2 = abs(xp - euclid(x_hat_p));

error = [error1(:) error2(:)];