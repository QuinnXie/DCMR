function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MRMinimax)
    # if prob.solution.solnError.errorTNIR[end] >= meshRefi.resiThreshold
    if (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold)
        numMesh = 0;
        N = alg.numInterval
        # First meshRefinement to find order reduction r
        errorMean = prob.solution.solnError.errorMean;
        numInsert = Int64.(ones(N))
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)
        solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;

        # order of accuracy
        p = typeof(alg) <: Trapezoidal ? 2 : 4;

        I₁ = 1; # number of points added in interval a
        r̂ = []
        for i in eachindex(numInsert)
            θ = errorMean[i]; # oldGridError
            a = i + i - 1; # index of new mesh
            η = sum(prob.solution.solnError.errorMean[a:(a+I₁)])/(1+I₁); # newGridError
            r̂ = [r̂; round(p + 1 - log(θ/η)/log(1+I₁))]
        end
        # @show r̂
        r = max.(0,min.(r̂,p)); 
        # r = Int64(round(sum(r)/N)) # order of reduction r for I₁ = 1
        # println("r = ", r)
    end

    # while (prob.solution.solnError.errorTNIR[end] >= meshRefi.resiThreshold)# && numMesh <=50)
    while (findmax(prob.solution.solnError.errorMean)[1] > meshRefi.resiThreshold)
        numInsert = meshRefinement(r, prob.solution.solnError.errorMean, meshRefi.resiThreshold, prob.solution.solnError.errorTime, alg, meshRefi);
        prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)

        solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
        numMesh = numMesh + 1;
        # @show "MRMinimaxError"
    end
    return prob.solution, alg
end

function meshRefinement(r::AbstractArray, errorMean::AbstractArray, resiThreshold::AbstractFloat, errorTime::AbstractArray, alg::Collocation, meshRefi::MRMinimax)
    N = alg.numInterval
    M = N + 1; # number of knots
    I = zeros(N); # store number of points to be added to each interval
    ΔM = 0; # count the total number of points added
    M₁ = 5; # maximum number of points added to a single interval
    κ = 1/10; # auxiliary constant 
    δ = resiThreshold; # limit of maximum error

    initial_time_segment = errorTime[end]/length(r);
    # order of accuracy
    # Trapezoidal - p = 2; HermiteSimpson - p = 4
    p = typeof(alg) <: Trapezoidal ? 2 : 4;

    ϵ = errorMean
    ϵₐ, a = findmax(ϵ)
    while ~((ΔM >= min(M₁,κ*M) && ϵₐ <= δ && I[a] == 0) || (ϵₐ <= κ*δ && 0 <= I[a] <= M₁) || (ΔM == M - 1) || (maximum(I)==M₁))
        # @show mod(errorTime[a],initial_time_segment)
        # @show Int64(floor.(errorTime[a]/initial_time_segment)+1)
        rₑ = r[Int64(floor.(errorTime[a]/initial_time_segment)+1)]

        # add a point to interval α
        I[a] = I[a] + 1 # number of points added to each interval
        ΔM = ΔM + 1 # total number of points added

        # update the predicted error for interval α
        ϵ[a] = ϵₐ*(1/(1+I[a]))^(p-rₑ[1]+1)

        # Recompute the maximum error and corrseponding interval α
        ϵₐ, a = findmax(ϵ)
        # @show ϵₐ
        # @show a
    end
    # @show I
    # @show N
    return I
end

# Another way to find reduction order r
# function meshRefinement(prob::DOProblem, alg::Collocation, meshRefi::MRMinimax)
#     numMesh = 0;

#     # First meshRefinement to find order reduction r
#     errorMean = prob.solution.solnError.errorMean;
#     maxError = findmax(errorMean)
#     numInsert = errorMean .>= maxError[1]
#     prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)
#     solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
#     numMesh = numMesh + 1;

#     # order of accuracy
#     p = typeof(alg) <: Trapezoidal ? 2 : 4;

#     I₁ = 1; # number of points added in interval a
#     a = maxError[2]; # index of maximum error
#     θ = errorMean[a]; # oldGridError
#     η = sum(prob.solution.solnError.errorMean[a:(a+I₁)])/(1+I₁); # newGridError
#     r̂ = round(p + 1 - log(θ/η)/log(1+I₁));
#     r = Int64(max(0,min(r̂,p))); # order of reduction r for I₁ = 1

#     while (prob.solution.solnError.errorTNIR[end] >= meshRefi.resiThreshold)# && numMesh <=50)
#         numInsert = meshRefinement(r, prob.solution.solnError.errorMean, meshRefi.resiThreshold, alg, meshRefi);
#         prob.solution, alg = meshRefinement(numInsert, prob, alg, meshRefi)

#         solve(prob, alg, meshRefi); # println("numInterval = ", alg.numInterval)
#         numMesh = numMesh + 1;
#         # @show meshMaxError = findmax(prob.solution.solnError.errorMean);
#         # @show meshMinError = findmin(prob.solution.solnError.errorMean);
#     end
#     return prob.solution, alg
# end