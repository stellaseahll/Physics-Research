using QuantumOptics
using Arpack
import Base:setproperty! #must be imported to overload it later for the present objects
import SparseArrays: SparseMatrixCSC, diag, diagm

#global parameters 
verbose = true; #whether to print info for relevant computation steps
wp = 1.0; #pointer frequency: sets the frequency/time units for all the following!
Bs = SpinBasis(1//2); #qubit basis defined for all the following 
P_ee = sigmap(Bs)*sigmam(Bs); #projector on excited state 
P_gg = sigmam(Bs)*sigmap(Bs); #projector on ground state 

#type abbreviations
typeBp = FockBasis;
typeB = CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}};
typeDO = DenseOperator{CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},Array{Complex{Float64},2}};
typeSO = SparseOperator{CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},
            SparseMatrixCSC{Complex{Float64},Int64}};
typeSSO = SparseSuperOperator{Tuple{CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}}},
            Tuple{CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}},CompositeBasis{Tuple{SpinBasis{1//2},FockBasis}}}, 
            SparseMatrixCSC{Complex{Float64},Int64}}

#useful simple functions for everyone
#compute nbar from T and transition frequency w
getNfromT(w::Float64, T::Float64) = 1.0/(exp(w/T)-1.0)




##########################
#The qubit-pointer models#
##########################


#abstract qubit + oscillator pointer model: Assuming qubit + pointer + displacement coupling, as well as hot qubit bath and cold pointer bath
abstract type qubitPointerModel end

mutable struct drivenQubitPointerModel <: qubitPointerModel
    #Assumes coherent driving with a driving freq wd. Hamiltonian is then given in rotating frame (unique fields: wd, G)
    # We will use a correlated basis defined by the vectors |m,n> where |g,n> = |g> ⊗ D|n> and |e,n> = |e> ⊗ D' |n>
    #  with D the displacement operator to +g/2ω
    # This can be implemented by an abstract tensor basis, but then be careful with the definition of spin- or pointer-only operators!
    
    N::Int #max Fock number  
    hotGlobal::Bool #whether hot bath is global (default=true)
    coldGlobal::Bool #whether cold bath is global (default=true)
    eigTol::Float64 #tol for arpack eigs solver (default = 1e-9)
    eigIter::Int #max iterations for arpack eigs solver (default = 500)
    
    isCooked::Bool #whether Hamiltonian is computed and consistent with parameters (internal check)
    isBathed::Bool #whether hot & cold baths are computed and consistent with parameters (internal check)
    isFed::Bool #total Liouvillian is computed and consistent with parameters (internal check)
   
    #all frequencies in units of constant wp (see above)
    ws::Float64 #spin frequency 
    wd::Float64 #driving frequency (optimal: red-detune by g^2/wp )
    g::Float64 #qubit-pointer displacement coupling 
    G::Float64 #effective driving rate 
    kh::Float64 #hot bath rate
    kc::Float64 #cold bath rate
    Th::Float64
    Tc::Float64 #bath temperatures in units hbar*omega{spin,pointer}/kB

    #basic bases, operators and stuff
    Bp::typeBp
    B::typeB

    DD::DenseOperator{FockBasis,FockBasis,Array{Complex{Float64},2}} #displacement operator D(g/wp) ONLY on Fock Hilbert space, is dense 
    SP::typeSO #sigma_+ operator in the correlated basis representation
    
    #all properties below get initialized in separate function 
    H::typeSO # system Hamiltonian in wd-rotating frame
  
    rhoSS::typeDO #the steady state of the system, will be dense 

    Lh::typeSSO #hot bath Liouvillian
    Lc::typeSSO #cold bath Liouvillian
    U::typeSSO #unitary Liouvillian    
    L::typeSSO #total Liouvillian

    function drivenQubitPointerModel(nmax, ws0, wd0, g0, G0, kh0, kc0, Th0, Tc0; 
                hotGlobal::Bool=true, coldGlobal::Bool=true, eigTol::Float64=1e-9, eigIter::Int=500, shitHappens::Bool=true )
        if nmax <= 0 || !isa(nmax,Int)
            error("max Fock number must be a positive integer")
        end
        Bp = FockBasis(nmax)
        bb = tensor(Bs,Bp)
        #incomplete initialization of object. Operators and Liouvillians still undefined
        m = new(nmax, hotGlobal, coldGlobal, eigTol, eigIter, false, false, false, 
                    ws0, wd0, g0, G0, kh0, kc0, Th0, Tc0, Bp, bb)
        #complete object by computing all operators and superoperators
        if shitHappens 
            cookHam!(m) #Hamiltonian and important auxiliary operators
            takeBath!(m) #cold and hot bath Liouvillians 
            feedLiouvillian!(m) #total Liouvillian and steady state 
        end
        return m 
    end

end

mutable struct randMeasQubitPointerModel <: qubitPointerModel
    #Assumes random projective pointer measurement at rate gamma + conditional feedback if LEFT. 
    #  LEFT means projection to |e,n> displ. Fock states with energy lower than potential at x=0.
    # We will use a correlated basis defined by the vectors |m,n> where |g,n> = |g> ⊗ D|n> and |e,n> = |e> ⊗ D' |n>
    #  with D the displacement operator to +g/2ω
    # This can be implemented by an abstract tensor basis, but then be careful with the definition of spin- or pointer-only operators!
    
    N::Int #max Fock number  
    hotGlobal::Bool #whether hot bath is global (default=true)
    coldGlobal::Bool #whether cold bath is global (default=true)
    eigTol::Float64 #tol for arpack eigs solver (default = 1e-9)
    eigIter::Int #max iterations for arpack eigs solver (default = 500)
    
    isCooked::Bool #whether Hamiltonian is computed and consistent with parameters (internal check)
    isBathed::Bool #whether hot & cold baths are computed and consistent with parameters (internal check)
    isFed::Bool #total Liouvillian is computed and consistent with parameters (internal check)
   
    #all frequencies in units of constant wp (see above)
    ws::Float64 #spin frequency 
    g::Float64 #qubit-pointer displacement coupling 
    gamma::Float64 #measurement rate  
    kh::Float64 #hot bath rate
    kc::Float64 #cold bath rate
    Th::Float64
    Tc::Float64 #bath temperatures in units hbar*omega{spin,pointer}/kB

    #basic bases, operators and stuff
    Bp::typeBp
    B::typeB

    DD::DenseOperator{FockBasis,FockBasis,Array{Complex{Float64},2}} #displacement operator D(g/wp) ONLY on Fock Hilbert space, is dense 
    SP::typeSO #sigma_+ operator in the correlated basis representation
    
    #all properties below get initialized in separate function 
    H::typeSO # system Hamiltonian in wd-rotating frame
  
    rhoSS::typeDO #the steady state of the system, will be dense 

    Lh::typeSSO #hot bath Liouvillian
    Lc::typeSSO #cold bath Liouvillian
    Lmfb::typeSSO #meas+fb Liouvillian (LEFT)
    Lm::typeSSO #meas Liouvillian (RIGHT)
    U::typeSSO #unitary Liouvillian    
    L::typeSSO #total Liouvillian

    function randMeasQubitPointerModel(nmax, ws0, g0, gamma0, kh0, kc0, Th0, Tc0; 
                hotGlobal::Bool=true, coldGlobal::Bool=true, eigTol::Float64=1e-9, eigIter::Int=500, shitHappens::Bool=true )
        if nmax <= 0 || !isa(nmax,Int)
            error("max Fock number must be a positive integer")
        end
        Bp = FockBasis(nmax)
        bb = tensor(Bs,Bp)
        #incomplete initialization of object. Operators and Liouvillians still undefined
        m = new(nmax, hotGlobal, coldGlobal, eigTol, eigIter, false, false, false, 
                    ws0, g0, gamma0, kh0, kc0, Th0, Tc0, Bp, bb)
        #complete object by computing all operators and superoperators
        if shitHappens 
            cookHam!(m) #Hamiltonian and important auxiliary operators
            takeBath!(m) #cold and hot bath Liouvillians 
            feedLiouvillian!(m) #total Liouvillian and steady state 
        end
        return m 
    end

end




function setproperty!(obj::qubitPointerModel,name::Symbol,x)
    #overload m.name = x to do detect parameter changes. 
    setfield!(obj,name,x) #assign the new obj.name = value. Error, if property doesn't exist
    if name == :N 
        obj.Bp = FockBasis(x)
        obj.B = tensor(Bs,obj.Bp)
        obj.isCooked = false; obj.isBathed = false; obj.isFed = false;
        verbose && println("New Fock basis! Please cookHam! -> takeBath! -> feedLiouvillian!")
    elseif name == :g || name == :ws
        obj.isCooked = false; obj.isBathed = false; obj.isFed = false;
        verbose && println("New system! Please cookHam! -> takeBath! -> feedLiouvillian!")
    elseif name == :G || name == :wd
        obj.isCooked = false; obj.isFed = false;
        verbose && println("New driving strength! Please cookHam! -> feedLiouvillian!")
    elseif name == :kh || name == :kc || name == :Th || name == :Tc || name == :hotGlobal
        obj.isBathed = false; obj.isFed = false;
        verbose && println("New baths! Please takeBath! -> feedLiouvillian!")
    elseif name == :eigTol || name == :eigIter
        verbose && println("New shitting routine! Please feedLiouvillian! or shitState!")
    end
end




#############################################################
#important functions for the relevant results & computations#
#############################################################

function cookHam!(m::qubitPointerModel)
    #computes Hamiltonian and all necessary operators, then also unitary superoperator
    verbose && println("I love the smell of cooking ham...")
    m.DD = displace(m.Bp,m.g/wp) #Fock representation of displacement operator between left and right
    #We have a weird basis of conditionally displaced states, so sigma_+ and sigma_- must include relative displacements
    m.SP = tensor(sigmap(Bs), m.DD) #sigma_plus
    # m.H = (m.ws-m.wd)/2*tensor(sigmaz(Bs), identityoperator(m.Bp)) + wp*tensor(identityoperator(Bs), number(m.Bp)) + m.G*(m.SP + dagger(m.SP));
    m.H = Hamiltonian(m) #see definitions depending on model below
    m.U = spre(-1im*m.H) + spost(1im*m.H) #superoperator for unitary part
    m.isCooked = true
    return nothing 
end

#no driving, Schroedinger frame 
Hamiltonian(m::qubitPointerModel) = 
                m.ws/2*tensor(sigmaz(Bs), identityoperator(m.Bp)) + wp*tensor(identityoperator(Bs), number(m.Bp));
#coherent driving term, in rotating frame with driving frequency wd                
Hamiltonian(m::drivenQubitPointerModel) = 
                (m.ws-m.wd)/2*tensor(sigmaz(Bs), identityoperator(m.Bp)) + wp*tensor(identityoperator(Bs), number(m.Bp)) + m.G*(m.SP + dagger(m.SP));




function takeBath!(m::drivenQubitPointerModel)
    #computes cold and hot bath Liouvillians (hot local or (semi-)global, depending on hotGlobal variable)
    verbose && println("Taking a cold and hot bath...")
    if !m.isCooked
        println("Cannot take hot bath! Please first cookHam!")
        return nothing
    end
    coldBath!(m)
    hotBath!(m)
    m.isBathed = true
    return nothing     
end

function takeBath!(m::randMeasQubitPointerModel)
    #computes cold and hot bath Liouvillians (hot local or (semi-)global, depending on hotGlobal variable)
    verbose && println("Taking cold, hot and random measurement bath...")
    if !m.isCooked
        println("Cannot take hot and measurement bath! Please first cookHam!")
        return nothing
    end
    coldBath!(m)
    hotBath!(m)
    measBath!(m)
    m.isBathed = true
    return nothing     
end


function hotBath!(m::qubitPointerModel)
    if m.hotGlobal
        #CAUTION: global dissipator does not include driving Hamiltonian, so not strictly global for driven scenario,  
        # i.e. only valid for weak driving!
        m.Lh = SparseSuperOperator((m.B,m.B),(m.B,m.B))
        for n=(-m.N):m.N
            nh = getNfromT(m.ws+n*wp,m.Th)
            LOp = tensor(sigmap(Bs), SparseOperator( m.Bp, m.Bp, diagm(-n => diag(m.DD.data,-n)) ) )
            m.Lh += liouvillian(SparseOperator(m.B), [ dagger(LOp), LOp ]; rates=[ m.kh*(nh+1), m.kh*nh ] )
        end
    else #local dissipator using sigma_plus and sigma_minus as jump operators
        nh = getNfromT(m.ws,m.Th)
        m.Lh = liouvillian(SparseOperator(m.B), [ dagger(m.SP), m.SP ]; rates=[m.kh*(nh+1), m.kh*nh])
    end
    return nothing 
end

function coldBath!(m::qubitPointerModel)
    nc = getNfromT(wp,m.Tc)
    a = tensor(identityoperator(Bs), destroy(m.Bp)) #annihilation operator
    if !m.coldGlobal #local cold bath (unphysical for strong displacements)
        #undisplaced mode operator is displaced one - sigma_z-dependent displacement
        a -= tensor(m.g/2/wp * sigmaz(Bs),identityoperator(m.Bp))
    end
    m.Lc = liouvillian(SparseOperator(m.B), [ a, dagger(a) ]; rates=[m.kc*(nc+1), m.kc*nc])
    return nothing 
end

function measBath!(m::randMeasQubitPointerModel)
    #random projective measurement at rate gamma. If LEFT, then do bit flip
    # LEFT: For convenience, we define projector in left-displaced Fock basis of excited qubit. 
    #  We count in all levels below the value of the harmonic potential at x=0, assuming x=(a+a')/sqrt(2) i.e. Hp=wp*(P^2/2+X^2/2)
    # (naive proj. halfspace meas. would induce infinite heating in principle...)
    Nmax = Int( ceil( max( (m.g/wp)^2/4.0 - 0.5, 1 ) ) ) #cutoff for projector
    #projector is TIMES sqrt(gamma)
    Pp = SparseOperator(m.Bp,m.Bp,sparse(1:Nmax+1, 1:Nmax+1, sqrt(m.gamma), m.N+1,m.N+1)) #projector in Fock space 
    P = tensor(P_ee,Pp) + tensor( P_gg, dagger(m.DD) * ( Pp*m.DD ) )
    sxP = tensor(sigmap(Bs), Pp*m.DD)
    sxP = sxP + dagger(sxP) 
    m.Lmfb = liouvillian(SparseOperator(m.B), [sxP]) #Dissipator describing measuring LEFT and feedback
    m.Lm = liouvillian(SparseOperator(m.B), [P]) #dissipator describing measuring RIGHT w/o feedback 
    return nothing 
end




function feedLiouvillian!(m::qubitPointerModel)
    #computes total Liouvillian m.L and steady state m.rhoSS
    if m.isCooked && m.isBathed 
        verbose && println("Liouvillian want Ham? Wah so clean! I give you Ham!")
        m.L = Liouvillian(m) #see definitions depending on model below
        m.isFed = true
        shitState!(m) #computes steady state 
    else
        verbose && print("Liouvillian hungry ma? ")
        !m.isCooked && print("Paiseh, forgot to cookHam! ")
        !m.isBathed && print("Wah so smelly, need to takeBath! ")
        println()
    end
    return nothing     
end

#default total Liouvillian: unitary + hot + cold bath 
Liouvillian(m::drivenQubitPointerModel) = m.U + m.Lh + m.Lc;
Liouvillian(m::randMeasQubitPointerModel) = m.U + m.Lh + m.Lc + m.Lmfb + m.Lm;



function shitState!(m::qubitPointerModel)
    if m.isFed
        verbose && println("Alamak! Liouvillian shitting a big stinky steady state on our carpet!")
        m.rhoSS = steadystate.eigenvector(m.L; tol=m.eigTol, which = :SM, maxiter=m.eigIter) #uses arpack eigs solver. SM= find eigs of smallest magnitude
    else
        println("Cannot shitState! First must feedLiouvillian!")
    end
    return nothing 
end




function getWorkHeat(m::drivenQubitPointerModel)
    #returns three values W, Qh, Qc given the system's steady state
    #For heat need to add free qubit term
    #drivenQubitPointerModel classes should have driving frequency wd defined
    H_Schr = m.H + m.wd/2*tensor(sigmaz(Bs),identityoperator(m.Bp))
    W = imag( 2*m.G*m.wd * tr( tensor(sigmap(Bs), m.DD)*m.rhoSS ) )
    Qh = real( tr( H_Schr*(m.Lh*m.rhoSS) ) )
    Qc = real( tr( H_Schr*(m.Lc*m.rhoSS) ) )
    verbose && println("W+Qh+Qc = ",W+Qh+Qc)
    return W,Qh,Qc
end

function getWorkHeat(m::randMeasQubitPointerModel)
    #returns 4 values W, Qh, Qc, Qfb given the system's steady state
    #W is not easy to define here! We use the meas+fb term Lmfb and subtract the contribution Qfb from Lm due to pure projective measurement
    Wmfb = real( tr( m.H*(m.Lmfb*m.rhoSS) ) )
    Qfb = real( tr( m.H*(m.Lm*m.rhoSS) ) )
    Qh = real( tr( m.H*(m.Lh*m.rhoSS) ) )
    Qc = real( tr( m.H*(m.Lc*m.rhoSS) ) )
    verbose && println("W+Qh+Qc+2*Qfb = ",Wmfb+Qfb+Qh+Qc)
    return Wmfb-Qfb,Qh,Qc,Qfb
end


#PRECOMPILE THE SHIT 
precompile(drivenQubitPointerModel,(Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64) );
