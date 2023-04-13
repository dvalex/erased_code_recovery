function b = near(x,y, tol)
if nargin < 3, tol = 1.0e-8;end
b = abs(x-y) < tol;
end