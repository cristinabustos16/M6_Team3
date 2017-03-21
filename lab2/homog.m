function [ xh ] = homog( xe )
%Homog: homogeneous coordinates from euclidean
xh = [xe; ones(1, size(xe,2))]
end