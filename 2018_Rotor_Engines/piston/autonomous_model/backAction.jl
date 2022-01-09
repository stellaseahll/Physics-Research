#JULIA version of pendulumEngine2ModeDrive
function myhist(y,x)
    y = sort(y)
    miny = y[1];
    maxy = y[end];
    dy = (maxy-miny)/(x*1.0)
    edges = zeros(x)
    xo = zeros(x)
    for i = 1:x
        edges[i] = miny + i*dy;
        xo[i] = miny - 0.5*dy + i*dy;
    end
    # edges = linspace(miny,maxy,x+1);
    idx = 1;
    pdf = zeros(x+1)
    pdf[end] = length(y)
    for k = 1:length(y)
        while (y[k]>edges[idx] && idx<x)
            pdf[idx+1]= k-1
            idx += 1
        end
    end
    no = (pdf[2:end]-pdf[1:end-1])/length(y)
    return (no,xo)
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
  return nothing
end

function thermalState!(na, nth)
  #draws random numbers according to thermal distribution (classical)
  #if only one array, then just intensity na
  @inbounds @simd for i=1:length(na)
    na[i] = nth*randexp()
  end
  return nothing
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


function runPendulumEngineFast!(avgs, snaps, Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)

  dt = tmax/Nt

  NSave = Int(floor(Nt/nSave)) + 1
  NSnap = Int(floor(Nt/nSnap)) + 1

  #avgs must be fillable (NSave,m) array with m the number of outputs
  # m=14 outputs are (in this column order): t,x,x2,p,p2,na,na2,PW,PH
  if ~( size(avgs)[1]==NSave && size(avgs)[2]>=14 )
    println(size(avgs))
    error("avgs array provided doesnt have correct size!")
  end

  if NSnap>0
    #snaps array must be (Ntrials,NSnap*2)
    if ~( size(snaps)[1]==Ntrials && size(snaps)[2]>=1 )
      error("snaps array provided doesnt have correct size!")
    end
  end

  #running variables
  xt = zeros(Ntrials) #spin vector x
  yt = zeros(Ntrials) #spin vector y
  zt = -ones(Ntrials) #spin vector z
  Xt = zeros(Ntrials) #position
  Pt = zeros(Ntrials) #momentum
  # no need for extra arrays, just allocate floats
  na = 0.0
  # initial states
  initialRotorState!(Xt,Pt,k,mu,L0)
 

  @inbounds @simd for i=1:Ntrials
    #save initial averages
    sinX = sin(Xt[i])
    cosX = cos(Xt[i])      
    fH = (1.0+sinX)/2.0
    fC = (1.0-sinX)/2.0
    fH2 = fH*fH
    fC2 = fC*fC
    na = 0.5*(zt[i]+1)
    avgs[1,2] += Xt[i]
    avgs[1,3] += Xt[i]*Xt[i]
    avgs[1,4] += Pt[i]
    avgs[1,5] += Pt[i]*Pt[i]
    avgs[1,6] += na
    avgs[1,7] += xt[i]
    avgs[1,8] += yt[i]
    avgs[1,9] += zt[i]
    avgs[1,10] += g*na*sinX*Pt[i]/I
    avgs[1,11] += kappa*g*cosX*fH2*(nH-na)
    avgs[1,12] += kappa*g*cosX*fC2*(nC-na)
    avgs[1,13] += kappa*fH2*(nH-na)
    avgs[1,14] += kappa*fC2*(nC-na)
  end
  avgs[1,1] = 0.0 #TIME

  # #save only nSave-th step
  iSave = 1 #running index to check
  jSave = 2 #index in result array

  iSnap = 1

  gdt = g*dt
  kdt = kappa*dt
  kHdt = kappa/2.0*(2*nH+1)*dt
  kCdt = kappa/2.0*(2*nC+1)*dt
  HdW = sqrt(2*kappa*nH*dt)
  CdW = sqrt(2*kappa*nC*dt)

  for nt = 1:Nt
    #Euler steps
    #check whether result is saved (do it in same loop)
    if iSave==nSave
      iSave=1
      @inbounds @simd for i=1:Ntrials
        randHr = randn()
        randHi = randn()
        randCr = randn()
        randCi = randn()
        sinX = sin(Xt[i])
        cosX = cos(Xt[i])
        fH = (1.0+sinX)/2.0
        fC = (1.0-sinX)/2.0
        fH2 = fH*fH
        fC2 = fC*fC
        keffdt = (kHdt*fH2+kCdt*fC2)
        dfH = cosX/2.0
        dfC = -dfH
        dX = Pt[i]/I*dt
        dP = g*sinX*0.5*(1+zt[i])*dt - HdW*0.5*dfH*(xt[i]*randHi+yt[i]*randHr)
                  - CdW*0.5*dfC*(xt[i]*randCi+yt[i]*randCr)
        dx = -gdt*cosX*yt[i] - keffdt*xt[i]+ HdW*fH*zt[i]*randHr + CdW*fC*zt[i]*randCr 
        dy = gdt*cosX*xt[i] - keffdt*yt[i] - HdW*fH*zt[i]*randHi - CdW*fC*zt[i]*randCi
        dz = -2*keffdt*zt[i] - (fH2+fC2)*kdt + HdW*fH*(-xt[i]*randHr+yt[i]*randHi) 
                  + CdW*fC*(-xt[i]*randCr+yt[i]*randCi) 
        Xt[i] += dX
        Pt[i] += dP
        xt[i] += dx
        yt[i] += dy
        zt[i] += dz
        sinX = sin(Xt[i])
        cosX = cos(Xt[i])
        fH = (1.0+sinX)/2.0
        fC = (1.0-sinX)/2.0
        fH2 = fH*fH
        fC2 = fC*fC
        na = 0.5*(zt[i]+1)
        avgs[jSave,2] += Xt[i]
        avgs[jSave,3] += Xt[i]*Xt[i]
        avgs[jSave,4] += Pt[i]
        avgs[jSave,5] += Pt[i]*Pt[i]
        avgs[jSave,6] += na
        avgs[jSave,7] += xt[i]
        avgs[jSave,8] += yt[i]
        avgs[jSave,9] += zt[i]
        avgs[jSave,10] += g*na*sinX*Pt[i]/I
        avgs[jSave,11] += kappa*g*cosX*fH2*(nH-na)
        avgs[jSave,12] += kappa*g*cosX*fC2*(nC-na)
        avgs[jSave,13] += kappa*fH2*(nH-na)
        avgs[jSave,14] += kappa*fC2*(nC-na)
      end
      avgs[jSave,1] = nt*dt
      jSave += 1
    else
      iSave += 1
      @inbounds @simd for i=1:Ntrials
        randHr = randn()
        randHi = randn()
        randCr = randn()
        randCi = randn()
        sinX = sin(Xt[i])
        cosX = cos(Xt[i])
        fH = (1.0+sinX)/2.0
        fC = (1.0-sinX)/2.0
        fH2 = fH*fH
        fC2 = fC*fC
        keffdt = (kHdt*fH2+kCdt*fC2)
        dfH = cosX/2.0
        dfC = -dfH
        dX = Pt[i]/I*dt
        dP = g*sinX*0.5*(1+zt[i])*dt - HdW*0.5*dfH*(xt[i]*randHi+yt[i]*randHr)
                  - CdW*0.5*dfC*(xt[i]*randCi+yt[i]*randCr)
        dx = -gdt*cosX*yt[i] - keffdt*xt[i]+ HdW*fH*zt[i]*randHr + CdW*fC*zt[i]*randCr 
        dy = gdt*cosX*xt[i] - keffdt*yt[i] - HdW*fH*zt[i]*randHi - CdW*fC*zt[i]*randCi
        dz = -2*keffdt*zt[i] - (fH2+fC2)*kdt + HdW*fH*(-xt[i]*randHr+yt[i]*randHi) 
                  + CdW*fC*(-xt[i]*randCr+yt[i]*randCi) 
        Xt[i] += dX
        Pt[i] += dP
        xt[i] += dx
        yt[i] += dy
        zt[i] += dz
      end
    end

    # #check whether snapshot is saved
    # if iSnap<=NSnap && mod(nt,nSnap)==1
    #   snaps[:,iSnap] = Pt
    #   iSnap += 1
    #   #snaps[:,2*NSnap + iSnap] = nat
    # end
  end

  @inbounds @simd for s=2:14
    @inbounds @simd for z=1:NSave
      avgs[z,s] /= Ntrials
    end
  end

  return nothing
end

#wrapper if new result arrays must be produced
function runPendulumEngineFast(Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)
  #allocate result arrays
  NSave = Int(floor(Nt/nSave)) + 1
  NSnap = Int(floor(Nt/nSnap)) + 1
  avgs = zeros(NSave,14)
  snaps = zeros(Ntrials,NSnap)
  runPendulumEngineFast!(avgs,snaps,Ntrials, tmax, Nt,
            I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
            nSave, nSnap)
  return (avgs,snaps)
end

function parPendulumEngineFast(Nwork, Ntrials, tmax, Nt,
          I, g, gamma, kT, kappa, nC, nH, k, mu, L0,
          nSave, nSnap)
  #parallelized execution of runPendulumEngineFast
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
  NSnap = Int(floor(Nt/nSnap)) + 1
  results = [Future(pid) for pid in procs()]
  #give em the jobs, and wait until all are finished
  @sync begin
    for (idx,pid) in w
      @async results[idx] =
            @spawnat pid runPendulumEngineFast(Ntrials, tmax, Nt,
                  I, g, gamma, kT, kappa,
                  nC, nH, k, mu, L0, nSave, nSnap)
    end
  end
  #collect their products
  avg = zeros(NSave,14)
  avgtmp = zeros(NSave,14)
  snaps = zeros(Nwork*Ntrials,NSnap)
  @inbounds @simd for idx = 1:Nwork
    println("test")
    (avgtmp,snaps[((idx-1)*Ntrials+1):(idx*Ntrials),:]) = fetch(results[idx])
    avg += avgtmp/Nwork
  end
  return (avg, snaps)

end

