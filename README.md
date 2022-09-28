# DCMR
Direct Collocation with Mesh Refinement Methods for Solving a dynamic optimization problem (DOP).
As shown in the main function, *Main_Cart_Pole.jl*, we take [cart pole problem][2] as an simple example.
Then, there are many inputs for collocation and meshrefi functions to set properties of them. They could be selected manually before the arguments.
```julia
julia> alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod);
julia> meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod);
```
After defining the methods, function *solveMeshRefi* is used to solve the DOP and initial and final solution would be stored in variable *prob_alg*. Finally, with function *plotResults*, all results could be plotted and shown in pictures.
```julia
julia> alg = prob_alg = solveMeshRefi(alg, meshRefi);
julia> plotResults(prob_alg);
```

Detailed description of these methods such as collocation and mesh refinement methods could be found in [Betts][1].

[1]:https://doi.org/10.1137/1.9780898718577
[2]:https://epubs.siam.org/doi/10.1137/16M1062569
