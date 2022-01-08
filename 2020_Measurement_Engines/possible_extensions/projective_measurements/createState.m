function rho = createState(r,mode)

% creates state and basis vector depending on input
% mode = 1 (polar coordinate), mode = 2 (cartesian)
% default = mode 1
if (nargin == 1) 
    mode = 1;
end

if (mode == 2)
    nx = r(1);
    ny = r(2);
    nz = r(3);
else
    ny = r(1)*sin(r(2))*sin(r(3));
    nx = r(1)*sin(r(2))*cos(r(3));
    nz = r(1)*cos(r(2));
end

rho = [ (1+nz), (nx-1i*ny); (nx+1i*ny) (1-nz)]/2;
