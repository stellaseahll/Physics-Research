Testing of standard thermalization ME of HO but using x and p
# We can do the sanity checks here
using QuantumOptics
using PyPlot
include("myHelpers.jl")

#################### Frequencies ####################
ω = 1 #ho frequency
k = 0.01 #cold bath coupling freq

#################### Temperatures and Occupation Numbers ####################
T = [0.5:0.5:5;] #hot bath temperature
data = zeros(7,10)
for it = 1:length(T)
    println(it)
    n = 1/(exp(ω/T[it])-1)

    #################### Harmonic Oscillator ####################
    # Define position basis
    # Note: for numerical stability, [xmin xmax] ~ [pmin pmax]
    Xmin = -10
    Xmax = 10
    N = 50
    dx = (Xmax-Xmin)/N
    bX = PositionBasis(Xmin,Xmax,N) #forms position basis
    bP = MomentumBasis(bX) #[pmin pmax] = [-pi/dx pi/dx]
    # Define Operators in momentum basis
    P = momentum(bX) # p op will be full in position basis
    X = position(bX) # x op will be sparse in position basis
    id = one(bX) # Identity for harmonic oscillator
    tfXP = transform(bX,bP) # Transform from p to x
    tfPX = transform(bP,bX) # Transform from x to p

    ####################  Hamiltonians ################

    H = ω/2 * (P*P+X*X)

    ####################  Dissipators and Rates ################
    a = (X + (1im * P))/sqrt(2) #a
    adg = (X - (1im * P))/sqrt(2) #adagger
    JumpOps = [a,adg]
    rateVec = [k*(n+1),k*n]

    ####################  Expectation values of operators ################
    function calc_expVal(t, ρ)
        p1 = real(expect(X, ρ)) #x value
        p2 = real(expect(X*X, ρ)) #x^2 value
        p3 = real(expect(P, ρ)) #p value
        p4 = real(expect(P*P, ρ)) #p^2 value
        p5 = real(expect(a*adg, ρ)) #N+1
        p6 = real(expect(adg*a, ρ)) #N
        return p1, p2, p3, p4, p5, p6
    end

    ####################  Initialisation ################

    X0 = 0.
    P0 = 0.2
    σx = 0.1 #initial width of wave packet
    ψin = gaussianstate(bX,X0,P0,σx) # Wave function in position bases
    ψin = normalize(ψin)

    ####################  Time Evolution of Coupling Hamiltonian ################
    # tf = 10/k
    # dt = tf/100
    # tVec = [0:dt:tf;]
    # tout, expVal = timeevolution.master(tVec, ψin, H, JumpOps; rates=rateVec, fout=calc_expVal)

    ####################  Process Output ################
    # data = zeros(7,length(tVec))
    # data[1,:] = tout
    # for i = 2:7
    #     data[i,:]= [p[i-1] for p=expVal]
    # end
    #
    # filename =  string("omega$(ω)_n$(n)_k$(k)_tf$(tf).csv")
    # saveMatrix(data,filename)

    ####################  Steady State ################
    tout, rhoSS = steadystate.master(H, JumpOps; rates=rateVec)
    p1,p2,p3,p4,p5,p6 = calc_expVal(0, rhoSS[2])
    data[:,it] = [p1,p2,p3,p4,p5,p6,n]


end

filename = string("steadySS_omega$(ω)_k$(k).csv")
saveMatrix(data,filename)
