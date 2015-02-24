function [X, Y, Z] = contoursToMatrix(x, C, y)
%[X, Y, Z] = contoursToMatrix(x, C [, y])
%
% Calculates a matrix of level values from a set of contour curves. More or less
% the inverse of creating contours from a matrix of surface height values.
%
% x - is a vector of X-axis coordinates at which the contours in C were sampled.
%
% C - is a matrix of contour values, one per row.
%     Note that the following must be true: length(C(n,:)) == length(x)
%     Contours are defined by the pair C(n,:) "lower values" and C(N+1-n,:)
%     "upper values", where N = size(C, 1);
%
% y - is a vector of Y-axis coordinates at which to sample the contours.
%
% X, Y - the coordinate values for each element in Z as produced by meshgrid.
%
% Z - a matrix of integer height values, essentially the index number of the
%     contour curve.

if ~ exist('y')
    y = linspace(min(min(C)), max(max(C)), length(x));
end
[X, Y] = meshgrid(x, y);
Z = zeros(size(X));
Nrows = size(C, 1);
for n = 1:Nrows/2
    c1 = C(n,:);
    c2 = C(Nrows+1-n,:);
    indices = find((Y >= repmat(c1, length(y), 1)) & (Y <= repmat(c2, length(y), 1)));
    Z(indices) = n;
end
