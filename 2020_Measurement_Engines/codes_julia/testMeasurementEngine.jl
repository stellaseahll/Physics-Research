Testing of Measurement Engine in x basis of HO
using QuantumOptics
# using PyPlot
# include("myHelpers.jl")

#################### Frequencies ####################
ωp = 0.1 #ho frequency
ωs = 1. #qubit frequency
g = 0.1 #coupling frequency
kh = 0.01 #hot bath coupling freq
kc = 0.01 #cold bath coupling freq
km = 0.01 #measurement rate

#################### Temperatures and Occupation Numbers ####################
Th = 1.0 #hot bath temperature
Tc = 0.01 #cold bath temperature
nh = 1/(exp(ωs/Th)-1)
nc = 1/(exp(ωp/Tc)-1)

#################### Harmonic Oscillator ####################
# Define position basis
k = 20 #factor away from equilibrium position ~\pm g/wp
Xmin = -20*g/ωp
Xmax = 20*g/ωp
N = 101
bX = PositionBasis(Xmin,Xmax,N) #forms position basis
bP = MomentumBasis(bX)
ProjP = zeros(N,N) #Projector to +x
ProjM = zeros(N,N) #Projector to -x
for i = Int((N-1)/2+1):N
    ProjP[i,i] = 1
end
for i = 1:Int((N-1)/2)
    ProjM[i,i] = 1
end
ProjP = DenseOperator(bX,ProjP)
ProjM = DenseOperator(bX,ProjM)

# Define Operators in momentum basis
P = momentum(bX) # p op will be full in position basis
X = position(bX) # x op will be sparse in position basis
idHO = one(bX) # Identity for harmonic oscillator

####################  Spin  ####################
bs = SpinBasis(1//2)
sx = sigmax(bs)
sy = sigmay(bs)
sz = sigmaz(bs)
sp = sigmap(bs)
sm = sigmam(bs)
idQ = one(bs) # Identity for qubit

####################  Overall Basis + Hamiltonian  ################
bXTot = bX ⊗ bs #overall X basis
bPTot = bP ⊗ bs #overall P basis
tfXP = transform(bX,bP) # Transform from p to x
tfPX = transform(bP,bX) # Transform from x to p
Hq = ωs/2 * (idHO ⊗ sz) # Qubit Hamiltonian
Xeff = (X ⊗ idQ) + g/ωp* (idHO ⊗ sz) #Effective X operator shifted due to qubit
Hho = ωp/2 * (P*P ⊗ idQ) + ωp/2 * (Xeff*Xeff)
Htot = Hho + Hq

####################  Dissipators and Rates ################
hotDown = idHO ⊗ sm
hotUp = idHO ⊗ sp
coldDown = (Xeff + (1im * P) ⊗ idQ)/sqrt(2) #b
coldUp = (Xeff - (1im * P) ⊗ idQ)/sqrt(2) #bdagger
measureP = ProjP ⊗ idQ #measurement feedback, measure +x do nothing on spin
measureM = ProjM ⊗ sx #measurement feedback, measure -x flip spin
JumpOps = [hotDown,hotUp,coldDown,coldUp,measureP,measureM]
rateVec = [kh*(nh+1),kh*nh,kc*(nc+1),kc*nc,km,km]

####################  Expectation values of operators ################
function calc_expVal(t, ρ)
    p1 = real(expect(idHO ⊗ sz, ρ)) #mean spin value
    p2 = real(expect(X ⊗ idQ, ρ)) #x value
    p3 = real(expect(X*X ⊗ idQ, ρ)) #x^2 value
    p4 = real(expect(P ⊗ idQ, ρ)) #p value
    p5 = real(expect(P*P ⊗ idQ, ρ)) #p^2 value
    p6 = real(expect(Xeff, ρ)) #shifted x value
    p7 = real(expect(Xeff*Xeff, ρ)) #shifted x^2 value
    p8 = real(expect(ProjP ⊗ idQ, ρ)) #projection on +x
    p9 = real(expect(ProjM ⊗ idQ, ρ)) #projection on -x, just for sanity check
    return p1, p2, p3, p4, p5, p6, p7, p8, p9
end

####################  Initialisation ################

X0 = 0.
P0 = 0.
σx = 0.01*g/ωp #initial width of wave packet
ψho = gaussianstate(bX,X0,P0,σx) # Wave function in position bases
ψs = spindown(bs) # Qubit in ground state
ψin = ψho ⊗ ψs # initial state

####################  Time Evolution of Coupling Hamiltonian ################
dt = 0.01
tf = 10
tVec = [0:dt:tf;]
tout, expVal = timeevolution.master(tVec, ψin, H, J; rates=rateVec, fout=calc_expVal)


####################  Process Output ################
# Reshape expectation values
exp_sz = [p[1] for p=expVal]
exp_X = [p[2] for p=expVal]
exp_X2 = [p[3] for p=expVal]
exp_P = [p[4] for p=expVal]
exp_P2 = [p[5] for p=expVal]
exp_Xeff = [p[6] for p=expVal]
exp_Xeff2 = [p[7] for p=expVal]
exp_ProjP = [p[8] for p=expVal]
exp_ProjM = [p[9] for p=expVal]
