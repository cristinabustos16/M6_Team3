function Hfit = compute_Hfit(corners)
    % Compute the affinity that is used to center the transformed image.
    % Every object that is going to be plotted in the transformed image,
    % must be transformed by this affinity previously.

    xmin = corners(1);
    xmax = corners(2);
    ymin = corners(3);
    ymax = corners(4);

    sistema_fit = [xmin, 0   , 1, 0;
               0   , ymin, 0, 1;
               xmax, 0   , 1, 0;
               0   , ymax, 0, 1];
    term_indep_fit = [0; 0; xmax - xmin + 1; ymax - ymin + 1];
    sol_fit = sistema_fit \ term_indep_fit;
    Hfit = [sol_fit(1), 0         , sol_fit(3);
            0         , sol_fit(2), sol_fit(4);
            0         , 0         , 1];
    
    return
end