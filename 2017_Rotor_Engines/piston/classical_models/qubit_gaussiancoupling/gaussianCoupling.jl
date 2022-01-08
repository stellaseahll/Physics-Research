#


#JULIA version of pendulumEngine2ModeDrive
function thermalState(N, nth)
  #draws random numbers according to thermal distribution (classical)
  phi0 = zeros(N)
  n0 = zeros(N)
  rand!(phi0)
  randexp!(n0)
  ar = sqrt(nth*n0).*cos(2*pi*phi0)
  ai = sqrt(nth*n0).*sin(2*pi*phi0)
  return (ar, ai)
end

function thermalState!(ar, ai, nth)
  #draws random numbers according to thermal distribution (classical)
  #arrays ar, ai to be filled must have same length
  @inbounds @simd for i=1:length(ar)
    phi0 = 2*pi*rand()
    n0 = nth*randexp()
    ar[i] = sqrt(n0)*cos(phi0)
    ai[i] = sqrt(n0)*sin(phi0)
  end
  #return (ar, ai)
end


function initialRotorState(N,k,x0,p0)
  #initial rotor state as in Rotor Heat Engine
  if k==0
    x = 2*pi*rand(N)
    p = p0*ones(N)
  elseif k<0
    x = x0*ones(N)
    p = p0*ones(N)
  else
    x = (1/sqrt(2*k))*randn(N) + x0
    p = sqrt(k/2)*randn(N) + p0
  end
  return (x, p)
end

function initialRotorState!(x,p,k,x0,p0)
  #initial rotor state as in Rotor Heat Engine
  @inbounds @simd for i=1:length(x)
    x[i] = x0
    p[i] = p0
    if k==0
      x[i] = 2*pi*rand()
    elseif k>0
      x[i] += (1/sqrt(2*k))*randn()
      p[i] += sqrt(k/2)*randn()
    end
  end
  #return (x, p)
end

function fHotCold!(fH,fC,angle)
   factor = 7.81067137848756498641
   k = 2.40206456421390335265 # max height = 2, area of fH2 = 3/8*2pi
   fH[1] = exp(k*cos(angle-pi/2.0))/factor
   fC[1] = exp(k*cos(angle+pi/2.0))/factor
end

function checkInputArrays(avgs,snaps,Nt,nSave,Ntrials,nSnap)

  NSave = Int(floor(Nt/nSave)) + 1
  NSnap = length(nSnap)

  #avgs must be fillable (NSave,m) array with m the number of outputs
  # m=13 outputs are (in this column order): t,x,x2,p,p2,ar,ai,na,PW,QH1, QC1, QH2, QC2
  if ~( size(avgs)[1]==NSave && size(avgs)[2]==13 )
    println(size(avgs))
    error("avgs array provided doesnt have correct size!")
  end

  if NSnap>0
    #snaps array must be (Ntrials,NSnap*6)
    if ~( size(snaps)[1]==Ntrials && size(snaps)[2]==6 )
      error("snaps array provided doesnt have correct size!")
    end
  end
  return (NSave,NSnap)
end

function runPendulumEngine2ModeDrive!(avgs, snaps, Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)

  dt = tmax/Nt

  (NSave,NSnap) = checkInputArrays(avgs,snaps,Nt,nSave,Ntrials,nSnap)

  #running variables
  nati = zeros(Ntrials)
  xt = zeros(Ntrials)
  pt = zeros(Ntrials)
  # no need for extra arrays, just allocate floats
  sinXti = 0.0
  cosXti = 0.0
  fHti = zeros(1)
  fCti = zeros(1)
  # initial states
  #thermalState!(art, ait, nC)
  #thermalState!(art, ait, nH)
  initialRotorState!(xt,pt,k,mu,L0)
  @inbounds @simd for i=1:Ntrials
    sinXti = sin(xt[i])
    cosXti = cos(xt[i])
    fHotCold!(fHti,fCti,xt[i])
    fH2ti = fHti[1]*fHti[1]
    fC2ti = fCti[1]*fCti[1]
    #save initial averages
    avgs[1,2] += xt[i]
    avgs[1,3] += xt[i]*xt[i]
    avgs[1,4] += pt[i]
    avgs[1,5] += pt[i]*pt[i]
    avgs[1,8] += nati[i]
    avgs[1,9] += g*nati[i]*sinXti*pt[i]/I
    avgs[1,10] += kappa*g*cosXti*fH2ti*(nH-nati[i])
    avgs[1,11] += kappa*g*cosXti*fC2ti*(nC-nati[i])
    avgs[1,12] += kappa*fH2ti*(nH-nati[i])
    avgs[1,13] += kappa*fC2ti*(nC-nati[i])
  end
  avgs[1,1] = 0.0 #TIME

  #save only nSave-th step
  iSave = 1 #running index to check
  jSave = 2 #index in result array

  iSnap = 1

  skHdW = sqrt(kappa*nH*dt/2)
  skCdW = sqrt(kappa*nC*dt/2)
  pH1 = nH*kappa*dt
  pC1 = nC*kappa*dt
  pH2 = (1+nH)*kappa*dt
  pC2 = (1+nC)*kappa*dt
  dampdW = sqrt(2*I*kT*gamma*dt)
  #dw = zeros(4)

  for nt = 1:Nt
    #Euler steps
    #check whether result is saved (do it in same loop)
    if iSave==nSave
      iSave=1
      @inbounds @simd for i=1:Ntrials
        # randn!(dw) NOT FASTER!
        sinXti = sin(xt[i])
        fHotCold!(fHti,fCti,xt[i])
        fH2ti = fHti[1]*fHti[1]
        fC2ti = fCti[1]*fCti[1]
        xt[i] += pt[i]/I*dt
        pt[i] += g*nati[i]*sinXti*dt
        if (nati[i]==0 && rand()<pH1*fH2ti+pC1*fC2ti)
          nati[i] = 1
          #r = rand()
          #if (r<0.25)
          #  pt[i] += 1
          #elseif (r<0.5)
          #  pt[i] -= 1
          #end
        elseif (nati[i]==1 && rand()<pH2*fH2ti+pC2*fC2ti)
          nati[i] = 0
          #r = rand()
          #if (r<0.25)
          #  pt[i] += 1
          #elseif (r<0.5)
          #  pt[i] -= 1
          #end
        end
        #save averages
        sinXti = sin(xt[i])
        cosXti = cos(xt[i])
        avgs[jSave,2] += xt[i]
        avgs[jSave,3] += xt[i]*xt[i]
        avgs[jSave,4] += pt[i]
        avgs[jSave,5] += pt[i]*pt[i]
        avgs[jSave,8] += nati[i]
        avgs[jSave,9] += g*nati[i]*sinXti*pt[i]/I
        avgs[jSave,10] += kappa*g*cosXti*fH2ti*(nH-nati[i])
        avgs[jSave,11] += kappa*g*cosXti*fC2ti*(nC-nati[i])
        avgs[jSave,12] += kappa*fH2ti*(nH-nati[i])
        avgs[jSave,13] += kappa*fC2ti*(nC-nati[i])
      end
      avgs[jSave,1] = nt*dt
      jSave += 1
    else
      iSave += 1
      @inbounds @simd for i=1:Ntrials
        # randn!(dw) NOT FASTER!
        sinXti = sin(xt[i])
        fHotCold!(fHti,fCti,xt[i])
        fH2ti = fHti[1]*fHti[1]
        fC2ti = fCti[1]*fCti[1]
        xt[i] += pt[i]/I*dt
        pt[i] += g*nati[i]*sinXti*dt
        if (nati[i]==0 && rand()<pH1*fH2ti+pC1*fC2ti)
          nati[i] = 1
          #r = rand()
          #if (r<0.25)
          #  pt[i] += 1
          #elseif (r<0.5)
          #  pt[i] -= 1
          #end
        elseif (nati[i]==1 && rand()<pH2*fH2ti+pC2*fC2ti)
          nati[i] = 0
          #r = rand()
          #if (r<0.25)
          #  pt[i] += 1
          #elseif (r<0.5)
          #  pt[i] -= 1
          #end
        end
      end
    end

    #check whether snapshot is saved
    if iSnap<=NSnap && nt==nSnap[iSnap]
      snaps[:,iSnap] = xt
      snaps[:,Nsnap + iSnap] = pt
      snaps[:,2*Nsnap + iSnap] = art
      snaps[:,3*Nsnap + iSnap] = ait
      snaps[:,4*Nsnap + iSnap] = art
      snaps[:,5*Nsnap + iSnap] = ait
    end
  end

  @inbounds @simd for s=2:13
    @inbounds @simd for z=1:NSave
      avgs[z,s] /= Ntrials
    end
  end

  return nothing
end

#wrapper if new result arrays must be produced
function runPendulumEngine2ModeDrive(Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)
  #allocate result arrays
  NSave = Int(floor(Nt/nSave)) + 1
  NSnap = length(nSnap)
  avgs = zeros(NSave,13)
  snaps = zeros(Ntrials,6*NSnap)
  runPendulumEngine2ModeDrive!(avgs,snaps,Ntrials, tmax, Nt,
            I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
            nSave, nSnap)
  return (avgs,snaps)
end



function parPendulumEngine2ModeDrive(Nwork, Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)
  #parallelized execution of runPendulumEngine2ModeDrive
  w = enumerate(procs()) #all workers ids with enumeration idx
  Nw = length(w)
  if Nw<Nwork
    println("ACHTUNg! Only ",Nw," angry French workers available! Reducing workload before they go on strike...")
    Nwork = Nw
  elseif Nw>Nwork
    println("ACHTUNg! Only ",Nwork," of ",Nw," workers needed. You could fire the others...")
  end
  print("-> ",Nwork," workers simulating ",Ntrials," trajectories each...")
  NSave = Int(floor(Nt/nSave)) + 1
  NSnap = length(nSnap)
  avgs = zeros(NSave,13)
  snaps = zeros(Nwork*Ntrials,6*NSnap)
  results = [Future(pid) for pid in procs()]
  @sync begin
    for (idx,pid) in w
      @async results[idx] =
            @spawnat pid runPendulumEngine2ModeDrive(Ntrials, tmax, Nt,
                  I, g, gamma, kT, kappa,
                  nC, nH, k, mu, L0, nSave, nSnap)
    end
  end
  avg = zeros(NSave,13)
  avgtmp = zeros(NSave,13)
  @inbounds @simd for idx = 1:Nwork
    (avgtmp,snaps[((idx-1)*Ntrials+1):(idx*Ntrials),:]) = fetch(results[idx])
    avg += avgtmp/Nwork
  end

  return (avg, snaps)

end