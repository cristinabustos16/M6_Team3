% The function "apply_H" that gets as input a homography and an image and 
% returns the image transformed by the homography.
function [transformedImage] = apply_H(I, H)

    % 1. compute the size of the transformed image

    sizeI = size(I);
        
    topLeft = double([0, 0, 1]');
    topRight = double([0, sizeI(2), 1]');
    bottomLeft = double([sizeI(1), 0, 1]');
    bottomRight = double([sizeI(1), sizeI(2), 1]');
    
    % apply H to corners: homogeneous coordinates
    cornersHomogeneous = [(H*topLeft)'; (H*topRight)'; (H*bottomLeft)'; (H*bottomRight)'];
    
    % corners: cartesian coordinates
    cornersCartesian = [cornersHomogeneous(:, 1)./cornersHomogeneous(:, 3) ...
                        cornersHomogeneous(:, 2)./cornersHomogeneous(:, 3)];

    % top left and bottom right new corners
    topLeft = min(cornersCartesian);
    bottomRight = max(cornersCartesian);


    % 2. initialize transformedImage

    transformedImage = uint8(zeros(round(bottomRight(1) - topLeft(1)), round(bottomRight(2) - topLeft(2)), 3));
        
    % 3. transformation of the pixels (backwards)

%     transposedI = permute(I,[2 1 3]);  % source: http://stackoverflow.com/questions/29707563/image-transpose-in-matlab

    H_inv = inv(H);

    for row=1:size(transformedImage, 1)
        for col=1:size(transformedImage, 2)
            p = [row col] + topLeft;
            p = double([p 1]');
            pH = H_inv * p;  % homogeneous coordinates
            pC = [round(pH(1)/pH(3)) round(pH(2)/pH(3))];  % cartesian coordinates
            
%             if  pC(2) > 0 && pC(1) > 0 && pC(2) <= size(transposedI,1) && pC(1) <= size(transposedI, 2)
            if  pC(2) > 0 && pC(1) > 0 && pC(1) <= size(I,1) && pC(2) <= size(I, 2)
%                 transformedImage(row, col, :) = transposedI(pC(2), pC(1), :);
                transformedImage(row, col, :) = I(pC(1), pC(2), :);
            end
        end
    end

end
