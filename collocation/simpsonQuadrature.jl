function simpsonQuadrature(fun, tLow::AbstractFloat, tUpp::AbstractFloat, n::Integer)
    nt = 2*n + 1;
    dt = (tUpp - tLow) / n;
    t = LinRange(tLow, tUpp, nt);

    nstate = size(fun(0),1)
    f = zeros(nt,nstate)
    for i = 1:nt
        f[i,:] .= abs.(fun(t[i]))
    end

    # Compute quadrature weights
    w = ones(nt,1)
    w[2:2:end] = 4 .* ones(length(w[2:2:end]))
    w[3:2:end-1] = 2 .* ones(length(w[3:2:end-1]))

    # Compute quadrature for each state (Absolute Local Error)
    if nstate == 1
        errorState = ((dt/6) * w' * f )[1]
    else
        errorState = ((dt/6) * w' * f )'
    end

    # # Compute quadrature for each interval
    # errorSqua = ((dt/6) * w' * sum(f.^2, dims = 2))[1]

    return errorState #, errorSqua
end