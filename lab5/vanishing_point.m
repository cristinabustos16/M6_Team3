% computes the vanishing point formed by the line that joins points xo1 and xf1
% and the line that joins points x02 and xf2
function [ v1 ] = vanishing_point(xo1, xf1, xo2, xf2)

    line1 = cross(xo1, xf1);
    line2 = cross(xo2, xf2);
    v1 = cross(line1, line2);

end