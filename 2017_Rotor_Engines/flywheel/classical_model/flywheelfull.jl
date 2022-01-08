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
  J = 0.5
  J2 = J*J
  #running variables
  # At = zeros(Ntrials) #spin vector x
  # Bt = zeros(Ntrials) #spin vector y
  zt1 = -J*ones(Ntrials) #spin vector z
  zt2 = -J*ones(Ntrials) #spin vector z
  xt1 = zeros(Ntrials) #spin vector x
  xt2 = zeros(Ntrials) #spin vector x
  yt1 = zeros(Ntrials) #spin vector y
  yt2 = zeros(Ntrials) #spin vector y
  
  Xt = zeros(Ntrials) #position
  Pt = zeros(Ntrials) #momentum

  # no need for extra arrays, just allocate floats
  na = 0.0
  # initial states
  initialRotorState!(Xt,Pt,k,mu,L0)
 

  @inbounds @simd for i=1:Ntrials
    #save initial averages  
    avgs[1,2] += Xt[i]
    avgs[1,3] += Xt[i]*Xt[i]
    avgs[1,4] += Pt[i]
    avgs[1,5] += Pt[i]*Pt[i]
    avgs[1,6] += zt1[i]
    avgs[1,7] += zt2[i]
  end
  avgs[1,1] = 0.0 #TIME

  # #save only nSave-th step
  iSave = 1 #running index to check
  jSave = 2 #index in result array

  iSnap = 1

  gdt = g/J*dt
  kdt = kappa/J*dt
  kncdt = 2*kappa/J*nC*dt
  knhdt = 2*kappa/J*nH*dt
  khcdt = kappa*(nH+nC)/J*dt
  hdW = sqrt(2*kappa*nH/J*dt)
  cdW = sqrt(2*kappa*nC/J*dt)

  for nt = 1:Nt
    #Euler steps
    #check whether result is saved (do it in same loop)
    if iSave==nSave
      iSave=1
      @inbounds @simd for i=1:Ntrials
        randH1 = randn()
        randH2 = randn()
        randC1 = randn()
        randC2 = randn()
        sinX = sin(Xt[i])
        cosX = cos(Xt[i])
        a1 = J2-zt1[i]*zt1[i]
        a2 = J2-zt2[i]*zt2[i]
        A = xt1[i]*xt2[i] + yt1[i]*yt2[i]
        B = xt1[i]*yt2[i] - yt1[i]*xt2[i]
        dX = Pt[i]/I*dt
        dP = gdt*sinX*A + gdt*cosX*B
        dz1 = -gdt*cosX*B - gdt*sinX*A - knhdt*zt1[i] -kdt*a1 + hdW*(randH1*yt1[i]-randH2*xt1[i])
        dz2 = gdt*cosX*B + gdt*sinX*A - kncdt*zt2[i]  -kdt*a2 + cdW*(randC1*yt2[i]-randC2*xt2[i])
        dx1 = gdt*cosX*zt1[i]*yt2[i] + gdt*sinX*zt1[i]*xt2[i] + hdW*randH2*zt1[i] + kdt*xt1[i]*zt1[i] - 0.5*knhdt*xt1[i]
        dx2 = gdt*cosX*zt2[i]*yt1[i] - gdt*sinX*zt2[i]*xt1[i] + cdW*randC2*zt2[i] + kdt*xt2[i]*zt2[i] - 0.5*kncdt*xt2[i]
        dy1 = -gdt*cosX*zt1[i]*xt2[i] + gdt*sinX*zt1[i]*yt2[i]- hdW*randH1*zt1[i] + kdt*yt1[i]*zt1[i] - 0.5*knhdt*yt1[i]
        dy2 = -gdt*cosX*zt2[i]*xt1[i] - gdt*sinX*zt2[i]*yt1[i] - cdW*randC1*zt2[i] + kdt*yt2[i]*zt2[i] - 0.5*kncdt*yt2[i]
        Xt[i] += dX
        Pt[i] += dP
        zt1[i] += dz1
        zt2[i] += dz2
        yt1[i] += dy1
        yt2[i] += dy2
        xt1[i] += dx1
        xt2[i] += dx2
        N1 = J/sqrt(zt1[i]*zt1[i] + yt1[i]*yt1[i] + xt1[i]*xt1[i])
        N2 = J/sqrt(zt2[i]*zt2[i] + yt2[i]*yt2[i] + xt2[i]*xt2[i])
        zt1[i] = zt1[i]*N1
        zt2[i] = zt2[i]*N2
        yt1[i] = yt1[i]*N1
        yt2[i] = yt2[i]*N2
        xt1[i] = xt1[i]*N1
        xt2[i] = xt2[i]*N2
        avgs[jSave,2] += Xt[i]
        avgs[jSave,3] += Xt[i]*Xt[i]
        avgs[jSave,4] += Pt[i]
        avgs[jSave,5] += Pt[i]*Pt[i]
        avgs[jSave,6] += zt1[i]
        avgs[jSave,7] += zt2[i]
      end
      avgs[jSave,1] = nt*dt
      jSave += 1
    else
      iSave += 1
      @inbounds @simd for i=1:Ntrials
        randH1 = randn()
        randH2 = randn()
        randC1 = randn()
        randC2 = randn()
        sinX = sin(Xt[i])
        cosX = cos(Xt[i])
        a1 = J2-zt1[i]*zt1[i]
        a2 = J2-zt2[i]*zt2[i]
        A = xt1[i]*xt2[i] + yt1[i]*yt2[i]
        B = xt1[i]*yt2[i] - yt1[i]*xt2[i]
        dX = Pt[i]/I*dt
        dP = gdt*sinX*A + gdt*cosX*B
        dz1 = -gdt*cosX*B - gdt*sinX*A - knhdt*zt1[i] -kdt*a1 + hdW*(randH1*yt1[i]-randH2*xt1[i])
        dz2 = gdt*cosX*B + gdt*sinX*A - kncdt*zt2[i]  -kdt*a2 + cdW*(randC1*yt2[i]-randC2*xt2[i])
        dx1 = gdt*cosX*zt1[i]*yt2[i] + gdt*sinX*zt1[i]*xt2[i] + hdW*randH2*zt1[i] + kdt*xt1[i]*zt1[i] - 0.5*knhdt*xt1[i]
        dx2 = gdt*cosX*zt2[i]*yt1[i] - gdt*sinX*zt2[i]*xt1[i] + cdW*randC2*zt2[i] + kdt*xt2[i]*zt2[i] - 0.5*kncdt*xt2[i]
        dy1 = -gdt*cosX*zt1[i]*xt2[i] + gdt*sinX*zt1[i]*yt2[i]- hdW*randH1*zt1[i] + kdt*yt1[i]*zt1[i] - 0.5*knhdt*yt1[i]
        dy2 = -gdt*cosX*zt2[i]*xt1[i] - gdt*sinX*zt2[i]*yt1[i] - cdW*randC1*zt2[i] + kdt*yt2[i]*zt2[i] - 0.5*kncdt*yt2[i]
        Xt[i] += dX
        Pt[i] += dP
        zt1[i] += dz1
        zt2[i] += dz2
        yt1[i] += dy1
        yt2[i] += dy2
        xt1[i] += dx1
        xt2[i] += dx2
        N1 = J/sqrt(zt1[i]*zt1[i] + yt1[i]*yt1[i] + xt1[i]*xt1[i])
        N2 = J/sqrt(zt2[i]*zt2[i] + yt2[i]*yt2[i] + xt2[i]*xt2[i])
        zt1[i] = zt1[i]*N1
        zt2[i] = zt2[i]*N2
        yt1[i] = yt1[i]*N1
        yt2[i] = yt2[i]*N2
        xt1[i] = xt1[i]*N1
        xt2[i] = xt2[i]*N2
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
    (avgtmp,snaps[((idx-1)*Ntrials+1):(idx*Ntrials),:]) = fetch(results[idx])
    avg += avgtmp/Nwork
  end
  return (avg, snaps)

end

