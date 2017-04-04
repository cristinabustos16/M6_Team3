function [ disparity_map ] = stereo_computation(left_image, right_image, min_disparity, max_disparity, window_size, matching_cost)
% Write a function called 'stereo_computation' that computes the disparity
% between a pair of rectified images using a local method based on a matching cost 
% between two local windows.

% - left image
% - right image
% - minimum disparity
% - maximum disparity
% - window size (e.g. a value of 3 indicates a 3x3 window)
% - matching cost (the user may able to choose between SSD and NCC costs)

[height,width] = size(left_image);%suppose left and right images are equally sized, and gray scale
disparity_map = zeros(size(left_image));
%window half is taken because we take a window around a pixel, border pixels 
%dont fulfill this restriction unless it is performed some kind of padding
if(mod(window_size,2)==0) %verify if it is even
    window_half_1 = window_size/2 - 1; 
    window_half_2 = window_size/2; 
else
    window_half_1 = floor(window_size/2); 
    window_half_2 = floor(window_size/2); 
end

start = ceil(window_size/2);
final = window_half_2;

for i = start : height - final
    for j = start : width - final
        cost_ssd = inf;
        cost_ncc = -inf;
        cost_bil = inf;
        
        block_left = left_image(i-window_half_1:i+window_half_2,j-window_half_1:j+window_half_2);
        
        for disp = min_disparity:max_disparity
            if(j+window_half_2+disp <= width)%height) %img boundaries constraints
               
                block_right = right_image(i-window_half_1 : i+window_half_2, ...
                                        j-window_half_1+disp : j+window_half_2+disp);%movement only in x-axis
                
                if strcmp(matching_cost, 'SSD')
                    current_cost = sum(sum((block_left-block_right).^2))/(window_size^2);
                    if current_cost < cost_ssd
                        cost_ssd = current_cost;
                        disparity_map(i,j) = disp;
                    end
                elseif strcmp(matching_cost, 'NCC')
                    I1 = sum(sum(block_left))./(window_size^2);
                    I2 = sum(sum(block_right))./(window_size^2);
                    sigmaI1 = sqrt(sum(sum((block_left - I1).^2)));
                    sigmaI2 = sqrt(sum(sum((block_right - I2).^2)));
                    current_cost = sum(sum( (block_left - I1).*(block_right - I2) )) / (sigmaI1*sigmaI2);
                    if current_cost > cost_ncc
                        cost_ncc = current_cost;
                        disparity_map(i,j) = disp;
                    end
                elseif strcmp(matching_cost, 'BILATERAL')
                    %paper values
                    gammac = 5;
                    gammap =17.5;
                    T = 10;
                    center = ceil(window_size/2);
                    I1 = block_left;
                    I2 = block_right;
                    I1_center = I1(center,center);
                    I2_center = I2(center,center);
                    num = 0;
                    den = 0;
                    for k = 1:window_size
                        for l = 1:window_size
                            dist = sqrt((k-center)^2 + (l-center)^2); 
                            wI1 = exp( - (abs(I1(k,l)-I1_center)/gammac) - (dist/gammap));
                            wI2 = exp( - (abs(I2(k,l)-I2_center)/gammac) - (dist/gammap));
                            c = min(I1_center - I2_center,T);
                            num = num + wI1*wI2*c;
                            den = den + wI1*wI2;
                        end
                    end
                    current_cost = num/den;
                    if current_cost < cost_bil
                        cost_bil = current_cost;
                        disparity_map(i,j) = disp;
                    end
                end
            end
        end
    end
end

end

