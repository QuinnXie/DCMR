function add_objective(model::Model, prob::DOProblem, alg::Trapezoidal{Continuous})
    N = alg.numInterval
    Δt = model.obj_dict[:Δt]
    u = model.obj_dict[:u]

    @NLobjective(model, Min, sum(Δt[j] / 2 * (u[j]^2 + u[j+1]^2) for j in 1:N))
end

function add_objective(model::Model, prob::DOProblem, alg::Trapezoidal{Discontinuous})
    N = alg.numInterval
    Δt = model.obj_dict[:Δt]
    u = model.obj_dict[:u]

    @NLobjective(model, Min, sum(Δt[j] / 2 * (u[j, 1]^2 + u[j, 2]^2) for j in 1:N))
end

function add_objective(model::Model, prob::DOProblem, alg::HermiteSimpson{Continuous})
    N = alg.numInterval
    Δt = model.obj_dict[:Δt]
    u = model.obj_dict[:u]

    @NLobjective(model, Min, sum(Δt[j] / 6 * (u[2*j-1]^2 + 4 * u[2*j]^2 + u[2*j+1]^2) for j in 1:N))
end

function add_objective(model::Model, prob::DOProblem, alg::HermiteSimpson{Discontinuous})
    N = alg.numInterval
    Δt = model.obj_dict[:Δt]
    u = model.obj_dict[:u]

    @NLobjective(model, Min, sum(Δt[j] / 6 * (u[j,1]^2 + 4 * u[j,2]^2 + u[j,3]^2) for j in 1:N))
end