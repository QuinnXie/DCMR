using PlotlyJS

open("plotResult\\medianTime.txt") do io
    median = readlines(io)
    global medianTime = parse.(Float64, median)
    println(medianTime)
end

## Median Running Time 
# start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
# MRMethodRun = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
# MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
collocationMethod = ["TRPuC"; "TRPuD"; "HSCuC"; "HSCuD"]

medianTimeCollocation = reshape(medianTime, 8, 4)
for j in 1:4
    # start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
    Term = ["defaultTerm", "earlyTerm"]
    collocationMethod = ["Trapezoidal", "Trapezoidal", "Hermite-Simpson", "Hermite-Simpson"]
    input_continuity = [Continuous(), Discontinuous(), Continuous(), Discontinuous()]
    MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]

    alg = Collocation(Ipopt.Optimizer(), 4, input_continuity[j], collocationMethod[j])
    meshRefi = MeshRefi(1e-6, warmStart(), lateTerm(), MRMethod[j])
    runningTime = medianTimeCollocation[:,j]
    pltRunTime = plotRunningTime(runningTime, Term, alg::Collocation, MRMethod)
end


# medianTimeMeshRefine = reshape(medianTime,4,16)
# for j in 1:4
#     start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
#     MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
#     collocationMethod = ["TRPuC"; "TRPuD"; "HSCuC"; "HSCuD"]
#     # MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
#     meshRefi = MeshRefi(1e-6, warmStart(), lateTerm(), MRMethod[j])
#     runningTime = medianTimeMeshRefine[1:4,j:4:16]
#     pltRunTime = plotRunningTime(runningTime, start, collocationMethod, meshRefi)
# end


Term = ["defaultTerm", "earlyTerm"]

function plotRunningTime(runningTime, Term, alg::Collocation, MRMethod)
    nTerm = length(Term)
    nMRMethod = length(MRMethod)

    if (typeof(alg) <: Trapezoidal{Continuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (TRPuC)"
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (TRPuD)"
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (HSCuC)"
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (HSCuD)"
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    pltRunTime = PlotlyJS.plot(
        [
            PlotlyJS.bar(name=Term[1], x=MRMethod, y=runningTime[1:nTerm:end-nTerm]),
            PlotlyJS.bar(name=Term[2], x=MRMethod, y=runningTime[2:nTerm:end-nTerm]),
            PlotlyJS.bar(name=Term[1], x=[MRMethod[end]], y=runningTime[end-nTerm+1:nTerm:end], yaxis="y2"),
            PlotlyJS.bar(name=Term[2], x=[MRMethod[end]], y=runningTime[end-nTerm+2:nTerm:end], yaxis="y2"),
        ],
        Layout(
            xaxis_title_text="Mesh Refinement Method",
            yaxis_title_text="Running Time (s)",
            yaxis2=attr(
                title="Running Time (s) (only for MRBisectionOne)",
                overlaying="y",
                side="right"
            ),
            title=attr(
                text=plot_title,
                y=0.95,
                x=0.5,
                xanchor= "center",
                yanchor= "top"
            ),
            legend=attr(
                x=1,
                y=-0.3,
                yanchor="bottom",
                xanchor="right",
                orientation="h",
                bordercolor="Black",
                borderwidth=2
            ),
        )
    )
    display(pltRunTime)
    # PlotlyJS.savefig(pltRunTime, format="png")
    return pltRunTime
end

# function plotRunningTime(runningTime, start, collocationMethod, meshRefi::MeshRefi)
    nTerm = length(Term)
    ncollocationMethod = length(collocationMethod)

    if typeof(meshRefi) <: MRMinimax
        MRType = "MRMinimax"
    elseif typeof(meshRefi) <: MREquidistributed
        MRType = "MREquidistributed"
    elseif typeof(meshRefi) <: MRBisectionMulti
        MRType = "MRBisectionMulti"
    elseif typeof(meshRefi) <: MRBisectionAll
        MRType = "MRBisectionAll"
    elseif typeof(meshRefi) <: MRBisectionOne
        MRType = "MRBisectionOne"
    else
        MRType = "$(typeof(meshRefi))"
    end

    plot_title = "Time for Solving DOP ($MRType)"

    pltRunTime = PlotlyJS.plot(
        [
            PlotlyJS.bar(name=Term[1], x=collocationMethod, y=runningTime[1:nTerm:end]),
            PlotlyJS.bar(name=Term[2], x=collocationMethod, y=runningTime[2:nTerm:end]),
        ],
        Layout(
            xaxis_title_text="Mesh Refinement Method",
            yaxis_title_text="Running Time (s)",
            title=attr(
                text=plot_title,
                y=0.95,
                x=0.5,
                xanchor= "center",
                yanchor= "top"
            ),
            legend=attr(
                x=1,
                y=-0.3,
                yanchor="bottom",
                xanchor="right",
                orientation="h",
                bordercolor="Black",
                borderwidth=2
            ),
        )
    )
    display(pltRunTime)
    # return pltRunTime
# end

