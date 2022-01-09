# Argument 

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--tstep"
            help = "tstep"
            arg_type = Float64
            default = 0.001
        "--nt"
            help = "nt"
            arg_type = Int
            default = 1000
        "--run"
            help = "run"
            arg_type = Int
            default = 0
        "--t"
            help = "t"
            arg_type = Float64
            default = 10.0
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
            default = 100.0
        "--g"
            help = "g"
            arg_type = Float64
            default = 10.0
        "--kappa"
            help = "kappa"
            arg_type = Float64
            default = 10.0
        "--kT"
            help = "kT"
            arg_type = Float64
            default = 0.0
        "--gamma"
            help = "gamma"
            arg_type = Float64
            default = 0.0
    end

    return parse_args(s)
end

parsed_args = parse_commandline()
println("Parsed args:")
for (arg,val) in parsed_args
    println("  $arg  =>  $val")
end


#SCRIPT that runs pendulumEngine simulation in parallel
# using PyPlot #for plotting

# How many cores (nW determines the total number of trials Nwork*Ntrials)
Nwork = 16;
#adds processes and loads function files to all of them
#Note that there is always one more process in the foreground, which remains idle/waits
addprocs(Nwork-1)
println("There are currently ",nprocs()," active processes")
# using ParallelDataTransfer
using DistributedArrays
@everywhere include("gaussianCoupling.jl")
include("myHelpers.jl")
#the parameters

Ntrials = parsed_args["nt"];
run_ = parsed_args["run"];
tmax = parsed_args["t"];
dt = parsed_args["tstep"];
Nt = Int(ceil(tmax/dt));
I_ = parsed_args["I"];
g_ = parsed_args["g"];
gamma_ = parsed_args["gamma"];
kT_ = parsed_args["kT"];
k = 10.0;
mu = pi/2;
L0 = 0.0;
kappa_ = parsed_args["kappa"];
nC = parsed_args["nc"];
nH = parsed_args["nh"];

nSave = 100;
nSnap = [];

#do a dry run with few trajectories to precompile
println("DRY RUN to precompile...")
@time parPendulumEngine2ModeDrive(Nwork, 100, 50, 1000,
          I_, g_, gamma_, kT_, kappa_, nC, nH, k, mu, L0,
          100, nSnap);


#now do the real one
println("FULL RUN...")
@time (avgs, snaps) = parPendulumEngine2ModeDrive(Nwork, Ntrials, tmax, Nt,
          I_, g_, gamma_, kT_, kappa_, nC, nH, k, mu, L0,
          nSave, nSnap)

#save averages
filename = string("DiffCouplingFn_Ntrials$(Ntrials)_tmax$(tmax)_I$(I_)_g$(g_)_kappa$(kappa_)_gamma$(gamma_)_kT$(kT_)_nH$(nH)_nC$(nC)_k$(k).csv")
saveMatrix(avgs, filename, "t, x, x2, p, p2, ar, ai, na, PW, QH1, QC1, QH2, QC2");

