# mesh refinement method type definition
abstract type MeshRefi end

# initial start type definition
abstract type MRStart end
struct coldStart <: MRStart end
struct onceWarmStart <: MRStart end
struct quasiWarmStart <: MRStart end
struct warmStart <: MRStart end

# termination type definition
abstract type MRTerm end
struct earlyTerm <: MRTerm end
struct lateTerm <: MRTerm end

## struct for mesh refinement method

# no mesh refinement
struct NoMR <: MeshRefi end

# mesh refinement Minimax error method
struct MRMinimax{S<:MRStart, T<:MRTerm}  <: MeshRefi 
    resiThreshold :: AbstractFloat
    start :: S
    termination :: T
end

# mesh refinement Bisection one method
struct MRBisectionOne{S<:MRStart, T<:MRTerm} <: MeshRefi 
    resiThreshold :: AbstractFloat
    start :: S
    termination :: T
end

# mesh refinement Bisection mutiple intervals method
struct MRBisectionMulti{S<:MRStart, T<:MRTerm} <: MeshRefi 
    resiThreshold :: AbstractFloat
    start :: S
    termination :: T
end

# mesh refinement Bisection all intervals method
struct MRBisectionAll{S<:MRStart, T<:MRTerm} <: MeshRefi 
    resiThreshold :: AbstractFloat
    start :: S
    termination :: T
end

# mesh refinement Bisection all intervals method
struct MREquidistributed{S<:MRStart, T<:MRTerm} <: MeshRefi 
    resiThreshold :: AbstractFloat
    start :: S
    termination :: T
end



function MeshRefi(
    resiThreshold::AbstractFloat = 1e-6,
    start::MRStart = warmStart(),
    termination::MRTerm = earlyTerm(),
    MRMethod::String = "MRMinimax"
)
    if MRMethod == "NoMR"
        return NoMR()
    elseif MRMethod == "MRBisectionOne"
        return MRBisectionOne(resiThreshold,start,termination)
    elseif MRMethod == "MRBisectionMulti"
        return  MRBisectionMulti(resiThreshold,start,termination)
    elseif MRMethod == "MRBisectionAll"
        return MRBisectionAll(resiThreshold,start,termination)
    elseif MRMethod == "MRMinimax"
        return  MRMinimax(resiThreshold,start,termination)  
    elseif MRMethod == "MREquidistributed"
        return  MREquidistributed(resiThreshold,start,termination)  
    else
        @error "Unexpected mesh refinement method '$(MRMethod)'"      
    end
end