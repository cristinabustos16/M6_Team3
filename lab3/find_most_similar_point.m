function x = find_most_similar_point(desc_ref, desc_N)
    % We use this function to locate SIFT points correpondent to the van in
    % the rest of the frames. Given a set of SIFT points, we look for the
    % one that is most similar to another SIFT point that is given.

    dist = zeros(1,size(desc_N,2));
    for i = 1:size(desc_N,2)
        dist(i) = norm(desc_N(:,i) - desc_ref);
    end
    
%     histogram(dist)

    [~, x] = min(dist);

return

end