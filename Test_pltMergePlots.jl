### merge Plots
pltTRPuC = plotResults(prob_alg_Group[1]);
pltTRPuD = plotResults(prob_alg_Group[5]);
pltHSSuC = plotResults(prob_alg_Group[9]);
pltHSSuD = plotResults(prob_alg_Group[13]);

# trajectory
pltTRPtraj = plot(
    pltTRPuC[3], 
    pltTRPuD[3],
    layout=grid(1,2),
    # linewidth=2,
    size=(800, 600),
    dpi = 1000
)
savefig(pltTRPtraj, "Cart_Pole_v4_plot\\plttraj\\pltTRPtraj")
savefig(pltTRPtraj, "Cart_Pole_v4_plot\\plttraj\\pltTRPtraj.svg")

pltHSStraj = plot(
    pltHSSuC[3], 
    pltHSSuD[3],
    layout=grid(1,2),
    # linewidth=2,
    size=(800, 600),
    dpi = 1000
)
savefig(pltHSStraj, "Cart_Pole_v4_plot\\plttraj\\pltHSStraj")
savefig(pltHSStraj, "Cart_Pole_v4_plot\\plttraj\\pltHSStraj.svg")

# time step
pltTimeStep = plot(
    pltTRPuC[2], 
    pltTRPuD[2],
    pltHSSuC[2], 
    pltHSSuD[2],
    layout=grid(2,2),
    size=(800, 600),
    dpi = 1000
)
savefig(pltTimeStep, "Cart_Pole_v4_plot\\pltStep\\pltTimeStep")
savefig(pltTimeStep, "Cart_Pole_v4_plot\\pltStep\\pltTimeStep.svg")

pltResτAll = plot(
    pltTRPuC[1], 
    pltTRPuD[1],
    pltHSSuC[1], 
    pltHSSuD[1],
    layout=grid(2,2),
    size=(850, 600),
    dpi = 1000
)
savefig(pltResτAll, "Cart_Pole_v4_plot\\pltResτ\\pltResτAll")
savefig(pltResτAll, "Cart_Pole_v4_plot\\pltResτ\\pltResτAll.svg")



touch("Cart_Pole_v4_plot\\runningTime\\warmStart.eps")
io  = open("warmStart.eps", "w")
    savefig(
        io,
        p::Plot;
        width::Union{Nothing,Int}=nothing,
        height::Union{Nothing,Int}=nothing,
        scale::Union{Nothing,Real}=nothing,
        format::String="png"
    )

