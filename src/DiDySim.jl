__precompile__(true)
module DiDySim
    using LightGraphs, SparseArrays, IterTools, StatsBase, GZip, GraphIO, Combinatorics, Formatting, Printf
    include("types.jl")
    #include("NetWiring.jl")
    include("NetWiring2.jl")
    #include("NetDynamics.jl")
    include("NetDynamics2.jl")
    #include("NetAnalysis.jl")
    include("NetAnalysis2.jl")
    export BNode,Net,BNode2,Net2,reg_function
    export init_node!, init_net!, rand_init_net!, evolve_node,evolve_node!, evolve_net!, obs_rand_cond!
    export derrida_map!, derrida_map2!,derrida_map3!, full_derrida_map!, derrida_map_mod!, derrida_map_mod2!, derrida_map_mod3!
    export knock_out!,knock_down!,fix_state!, get_basins!,find_attractor!,find_attractor2!,derrida_map_att3!
    export generate_RBN,generate_RBN2,load_net2, load_net, to_boolean!, read_metadata
end
