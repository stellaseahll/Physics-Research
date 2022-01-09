using QuantumOptics

using ArgParse
include("extractWork.jl")
include("myHelpers.jl")
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--t"
            help = "t"
            arg_type = Float64
            default = 500.0
        "--dt"
            help = "dt"
            arg_type = Float64
            default = 0.1
        "--nh"
            help = "nh"
            arg_type = Float64
            default = 1.0
        "--nc"
            help = "nc"
            arg_type = Float64
            default = 0.0
        "--I"
            help = "I"
            arg_type = Float64
            default = 10.0
        "--g"
            help = "g"
            arg_type = Float64
            default = 0.924
        "--kappa"
            help = "kappa"
            arg_type = Float64
            default = 1.0
        "--kT"
            help = "kT"
            arg_type = Float64
            default = 0.0
        "--gamma"
            help = "gamma"
            arg_type = Float64
            default = 0.0
        "--mmin"
            help = "mmin"
            arg_type = Int
            default = -50
        "--mmax"
            help = "mmax"
            arg_type = Int
            default = 160
    end
    return parse_args(s)
end

parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end

#----------------------QUBIT----------------------
qubitBasis = SpinBasis(1//2)
sm = sigmam(qubitBasis)
sp = sigmap(qubitBasis)
nh_ = parsed_args["nh"]; #bath occupation no.
nc_ = parsed_args["nc"];
ph = nh_/(2.0*nh_+1.0) #excited probability
pc = nc_/(2.0*nc_+1.0)
# Hh = 0.5*wh_*sigmaz(qubitBasis)
# Hc = 0.5*wc_*sigmaz(qubitBasis)
rhoh = ph*sp*sm + (1-ph)*sm*sp
rhoc = pc*sp*sm + (1-pc)*sm*sp
Iq = sp*sm+sm*sp
#----------------------ROTOR----------------------
m_min = parsed_args["mmin"];
m_max = parsed_args["mmax"];
m = m_min:m_max';
dimR = m_max - m_min +1
rotorBasis = GenericBasis(dimR)
k = 1
mu= 0.0
L0 = 0.0
r = zeros(1,dimR)
n = 0.0
MI = parsed_args["I"];
# Initialize matrix with_ zeros
lm=zeros(dimR,dimR) # exp(-i*phi)
lp=zeros(dimR,dimR) # exp(i*phi)
lz=zeros(dimR,dimR) # angular momentum
lzSquare=zeros(dimR,dimR) # angular momentum squared
rhor = zeros(dimR,dimR)
Ir = zeros(dimR,dimR)
# Compute matrix elements
for m = 1:dimR
    lz[m,m]=m_min+(m-1)
    Ir[m,m]=1
    lzSquare[m,m]=(m_min+m-1)*(m_min+m-1)
    if (m!=1)
        lm[m-1,m]=1 #Decrease in momentum
        lp[m,m-1]=1 #Increase in momentum
    end
    if (m==1-m_min)
        rhor[m,m]=1
    end

end
KEpartial = lzSquare/2.0/MI
cosPhi = 0.5*(lm+lp)
sinPhi = -0.5*1im*(lp-lm)
lm = sparse(DenseOperator(rotorBasis,lm))
lp = sparse(DenseOperator(rotorBasis,lp))
lz = sparse(DenseOperator(rotorBasis,lz))
lzSquare = sparse(DenseOperator(rotorBasis,lzSquare))
Ir = sparse(DenseOperator(rotorBasis,Ir))
cosPhi = sparse(DenseOperator(rotorBasis,cosPhi))
sinPhi = sparse(DenseOperator(rotorBasis,sinPhi))
rhor = DenseOperator(rotorBasis,rhor)
Hr = lzSquare/2.0/MI
KEfull = full(Iq⊗Hr⊗Iq)
KEfull = KEfull.data

#----------------------HAMILTONIANS----------------------
kappa_ = parsed_args["kappa"];
g_ = parsed_args["g"];

H0 = Iq⊗Hr⊗Iq
Hint = g_*(sm⊗lp⊗sp+sp⊗lm⊗sm)
H = H0 + Hint
H2 = (sp*sm)⊗Ir⊗Iq + Iq⊗Ir⊗(sp*sm)  
Htotal = full(H).data

#----------------------DISSIPATORS----------------------
a1 = sm⊗Ir⊗Iq
ad1 = sp⊗Ir⊗Iq
a3 = Iq⊗Ir⊗sm
ad3 = Iq⊗Ir⊗sp
#hot and cold baths
C1 = sqrt(kappa_*(nh_+1))*a1
C1d = sqrt(kappa_*nh_)*ad1
C2 = sqrt(kappa_*(nc_+1))*a3
C2d = sqrt(kappa_*nc_)*ad3

# #rotor friction
# kt = parsed_args["kT"];
# gamma_ = parsed_args["gamma"];
# gamma_ = 10^gamma_
# r1 = sqrt(2*kt*gamma_*MI)*(cosPhi - 1im*1.0/4.0/kt/MI * sinPhi *lz)
# r2 = sqrt(2*kt*gamma_*MI)*(sinPhi + 1im*1.0/4.0/kt/MI * cosPhi *lz)
# r1 = Iq⊗r1⊗Iq
# r2 = Iq⊗r2⊗Iq
# r1d = dagger(r1)
# r2d = dagger(r2)
J = [C1, C1d, C2, C2d]

#----------------------EXPECTATION VALUES----------------------
Qh = kappa_*(nh_+1)*(ad1*H2*a1 - 0.5*ad1*a1*H2 - 0.5*H2*ad1*a1) + kappa_*nh_*(a1*H2*ad1 - 0.5*a1*ad1*H2 - 0.5*H2*a1*ad1)
Qc = kappa_*(nc_+1)*(ad3*H2*a3 - 0.5*ad3*a3*H2 - 0.5*H2*ad3*a3) + kappa_*nc_*(a3*H2*ad3 - 0.5*a3*ad3*H2 - 0.5*H2*a3*ad3)
Qhload = kappa_*(nh_+1)*(ad1*H*a1 - 0.5*ad1*a1*H - 0.5*H*ad1*a1) + kappa_*nh_*(a1*H*ad1 - 0.5*a1*ad1*H - 0.5*H*a1*ad1)
Qcload = kappa_*(nc_+1)*(ad3*H*a3 - 0.5*ad3*a3*H - 0.5*H*ad3*a3) + kappa_*nc_*(a3*H*ad3 - 0.5*a3*ad3*H - 0.5*H*a3*ad3)
# W = (r1d*H*r1 - 0.5*r1d*r1*H - 0.5*H*r1d*r1) + (r2d*H*r2 - 0.5*r2d*r2*H - 0.5*H*r2d*r2)
ELz = Iq⊗lz⊗Iq
ELz2 = Iq⊗lzSquare⊗Iq

#---------------TIMESTEPS AND INITIALISATION----------------
t_ = parsed_args["t"]
dt_ = parsed_args["dt"]
timestep = [0:dt_:t_;]
rho0 = full(rhoh)⊗rhor⊗full(rhoc)

#---------------SOLVING MASTER EQUATION----------------
tout, rho = timeevolution.master(timestep,rho0,H,J)
#---------------EXPECTATION VALUES OF RUN----------------
data = zeros(17,length(timestep))
data[1,:] = timestep
data[5,:] = real(expect(Qh,rho))'
data[6,:] = real(expect(Qc,rho))'
data[8,:] = real(expect(Qhload,rho))'
data[9,:] = real(expect(Qcload,rho))'
data[10,:] = real(expect(ELz,rho))'
data[11,:] = real(expect(ELz2,rho))'
data[7,:] = data[10,:].^2/2.0/MI

#---------------EXTRACTABLE WORK----------------
for i = 1:length(timestep)
    rhoi = full(rho[i])
    rhoi_partial = ptrace(rho[i],[1,3]).data
    rhoi = rhoi.data
    data[2,i] = extractWork(rhoi,KEfull)
    data[3,i] = extractWork(rhoi,Htotal)
    data[4,i] = extractWork(rhoi_partial,KEpartial)  
    data[16,i] = real(rhoi[1,1])
    data[17,i] = real(rhoi[dimR,dimR])   
end

saveMatrix(data,"tmax$(t_)_dt$(dt_)_mmin$(m_min)_mmax$(m_max)_g$(g_)_I$(MI)_kappa$(kappa_)_nh_$(nh_)_nc_$(nc_).csv")
saveMatrix(full(rho[length(timestep)]),"State_tmax$(t_)_dt$(dt_)_mmin$(m_min)_mmax$(m_max)_g$(g_)_I$(MI)_kappa$(kappa_)_nh_$(nh_)_nc_$(nc_).csv")





