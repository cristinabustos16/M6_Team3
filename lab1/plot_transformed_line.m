% Compute and plot the transformed line:
function line_transformed = plot_transformed_line(line, I, H)

    % Translation applied for plotting:
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
    translation = min(cornersCartesian);
    
    % Transposed inverse of H:
    H_inv_tras = inv(H)';
    
    % Transformed line:
    line_transformed = H_inv_tras * line;
    
    % Plot the transformed line:
    t=1:0.1:1000;
    plot(t-translation(1), -(line_transformed(1)*t + line_transformed(3)) / line_transformed(2) - translation(2), 'r');

    return
    
end