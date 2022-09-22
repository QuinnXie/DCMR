function add_earlyTerm(model::Model, alg::Collocation, term::earlyTerm)
    N = alg.numInterval

    # with early termination
    # tol = max(1e-3/N, 1e-8)
    # tol = max(4*1e-3/N/N, 1e-8)
    # tol = max(16*1e-3/N^3, 1e-8)
    # tol = max(64*1e-3/N^4, 1e-8)

    tol = max(1e-2/N, 1e-8)
    # tol = max(4*1e-2/N/N, 1e-8)
    # tol = max(16*1e-2/N^3, 1e-8)
    # tol = max(64*1e-2/N^4, 1e-8)

    # @show tol
    # println("tol = ", tol)
    set_optimizer_attribute(model, "tol", tol)

    return tol;
end


function add_earlyTerm(model::Model, alg::Collocation, term::lateTerm)
    tol = 1e-8
    set_optimizer_attribute(model, "tol", tol);
end

