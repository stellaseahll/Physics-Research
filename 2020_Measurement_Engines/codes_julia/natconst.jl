module natconst

#just exporting the usefull natural constants in SI units as floats
export c, eV, h, hbar, eps0, my0, kB, amu, NA, me, G, a0, mbar, meV, Debye


const c = 299792458.0;
const eV = 1.60217733e-19;
const h = 6.6260755e-34;
const hbar = 1.05457266e-34;
const eps0 = 8.854187871e-12;
const my0 = 4*pi*1e-7;
const kB = 1.380658e-23;
const amu = 1.6605402e-27;
const NA = 6.0221367e23;
const me=9.1093897e-31;

const G = 6.67384e-11; #+/- 0.00080e-11 m^3/kg/s^2

const a0=4*pi*eps0*hbar^2/me/eV^2; #Bohr radius 5.29177249e-11; 

# andere Einheiten in SI
const mbar=100.0;
const barn=1.0e-28;
const meV=1.0e-3*eV;
const Debye=3.33564e-30;

mutable struct myStruct 
    x::Float64
    y::Float64
    z::Complex{x}
end

end #module

