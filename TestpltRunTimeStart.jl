using PlotlyJS

include("CollocationMeshRefine.jl")
include("cartPoleDynamics.jl")
include("plotResult\\plotRunTimeStart.jl")

open("plotResult\\medianTimeStart.txt") do io
    median = readlines(io)
    global medianTime = parse.(Float64, median)
    println(medianTime)
end

## Median Running Time 
# start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
# MRMethodRun = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
# MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
collocationMethod = ["TRPuC"; "TRPuD"; "HSSuC"; "HSSuD"]

medianTimeCollocation = reshape(medianTime, 16, 4)
for j in 1:4
    start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
    collocationMethod = ["Trapezoidal", "Trapezoidal", "Hermite-Simpson", "Hermite-Simpson"]
    input_continuity = [Continuous(), Discontinuous(), Continuous(), Discontinuous()]
    MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]

    alg = Collocation(Ipopt.Optimizer(), 4, input_continuity[j], collocationMethod[j])
    meshRefi = MeshRefi(1e-6, warmStart(), lateTerm(), MRMethod[j])
    runningTime = medianTimeCollocation[:,j]
    pltRunTime = plotRunningTime(runningTime, start, alg, MRMethod)
end


medianTimeMeshRefine = reshape(medianTime,4,16)
for j in 1:4
    start = ["coldStart"; "onceWarmStart"; "quasiWarmStart"; "warmStart"]
    MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
    collocationMethod = ["TRPuC"; "TRPuD"; "HSSuC"; "HSSuD"]
    # MRMethod = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
    meshRefi = MeshRefi(1e-6, warmStart(), lateTerm(), MRMethod[j])
    runningTime = medianTimeMeshRefine[1:4,j:4:16]
    pltRunTime = plotRunningTime(runningTime, start, collocationMethod, meshRefi)
    
end