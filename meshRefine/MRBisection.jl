function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MRBisectionAll)
    numMesh = 0;
    # while (prob.solution.solnError.errorTNIR[end]>=meshRefi.resiThreshold && numMesh <=20)
    while (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold && numMesh <=20)
        numInsert = Int64.(ones(alg.numInterval))
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)
            
        solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;
    end
    return prob.solution, alg
end

function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MRBisectionMulti)
    numMesh = 0;
    # while (prob.solution.solnError.errorTNIR[end]>=meshRefi.resiThreshold)
    while (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold)
        numInsert = prob.solution.solnError.errorMean.>meshRefi.resiThreshold
        # numInsert = prob.solution.solnError.errorMean .> (meshRefi.resiThreshold/alg.numInterval)
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)

        solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;
    end
    return prob.solution, alg
end

function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MRBisectionOne)
    numMesh = 0;
    # while (prob.solution.solnError.errorTNIR[end]>=meshRefi.resiThreshold)
    while (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold)
        maxError = findmax(prob.solution.solnError.errorMean)
        numInsert = prob.solution.solnError.errorMean .>= maxError[1]
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)
        solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;
    end
    return prob.solution, alg
end
