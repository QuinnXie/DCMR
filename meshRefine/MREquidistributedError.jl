function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MREquidistributed)
    numMesh = 0;
    while (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold)# && numMesh <=50)
        numInsert = meshRefinement(alg.numInterval, prob.solution.solnError.errorMean, meshRefi);
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi);

        soln, alg = solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;
    end
    return prob.solution, alg
end

function meshRefinement(N::Integer, errorMean::AbstractArray, meshRefi::MREquidistributed)
    M = N + 1 # number of knots

    E = cumsum([0; errorMean]) # monotonically increasing error
    m̂ = M + 1 # maximum number of points added totally is M - 1 (number of intervals)
    Ê = 0:sum(E[end])/(m̂-1):sum(E[end]) # error uniform distribution
    # @show E[end]
    # @show Vector(Ê)
    M₁ = 5 # the maximum number of points could be added for single interval 
    I = zeros(N) # number of points added in each interval

    # compare the uniform distribution error with the original error 
    # to find number of points added to each interval
    for i in 1:M-1
        for j in 2:m̂-1
            if  E[i] <= Ê[j] <= E[i+1]
                I[i] = I[i] + 1
            end
        end
        I[i] = min(I[i],M₁) # limit maximum points added to single interval
    end
    # @show sum(I)
    return I # number of points in each index
end