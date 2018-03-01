import GZip
using StatsBase
using IterTools
using LightGraphs
using Formatting
#derrida_map calculates derrida map for a net up to a  H hamming distance
function derrida_map!(netstate::Vector{BNode},nrand_init::Integer=prod([length(sigma.state_range) for sigma in netstate]);tau::Integer=1,Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netstate]))
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,netstate)
    N = length(index_nodes_to_flip)
    1<=Hmax<=N || error("Hmax must be less or equal than network size")
    Hmap = Tuple{Float64,Float64}[]
    pairs = [Dict{AbstractString,AbstractString}() for i in 1:Hmax]
    size_space_state = prod([length(sigma.state_range) for sigma in netstate])
    if nrand_init>size_space_state
        nrand_init = size_space_state
    end
    # tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    nredundant = 0
    Xa = Vector{Integer}(length(netstate))
    Xb = Vector{Integer}(length(netstate))
    i = 0
    while i < nrand_init
        Xa = rand_init_net!(netstate)
        Xa_tag = join(map(string,round(Int,Xa)))
        for H in 1:Hmax
            already_visitted = true
            while already_visitted
                Xb = deepcopy(Xa)
                for j in sample(index_nodes_to_flip, H, replace = false, ordered = true)
                    Xb[j] = rand(filter(x -> x!=Xa[j],netstate[j].state_range))
                end
                Xb_tag = join(map(string,round(Int,Xb)))
                if in((Xa_tag => Xb_tag),pairs[H]) #|| in((Xb_tag => Xa_tag),pairs[H])
                    nredundant += 1
                else
                    pairs[H][Xa_tag] = Xb_tag
                    already_visitted = false
                end
            end
    		for t in 1:tau
                for k in eachindex(netstate)
                    init_node!(netstate[k],Xa[k])
                end
    			Xa = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
    		end
            for t in 1:tau
                for k in eachindex(netstate)
                    init_node!(netstate[k],Xb[k])
                end
    			Xb = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
    		end
            push!(Hmap,(H/N,sum(Xa.!=Xb)/N))
            if mod(i,10000)==0
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
            end
        end
        i += 1
    end
    println("$nredundant revisitted configuration pairs were discarded")
    return reshape(Hmap,Hmax,nrand_init)
end

function derrida_map2!(netstate::Vector{BNode},nrand_init::Integer=prod([length(sigma.state_range) for sigma in netstate]);tau::Integer=1,Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netstate]))
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,netstate)
    N = length(index_nodes_to_flip)
    1<=Hmax<=N || error("Hmax must be less or equal than network size")
    Hmap = Tuple{Float64,Float64}[]
    pairs = [Dict{AbstractString,AbstractString}() for i in 1:Hmax]
    # tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    nredundant = 0
    Xa = Vector{Integer}(length(netstate))
    Xb = Vector{Integer}(length(netstate))
    Xatp1 = Vector{Integer}(length(netstate))
    Xbtp1 = Vector{Integer}(length(netstate))
    range_lengths = [length(sigma.state_range) for sigma in netstate]
    state_space_size = prod(range_lengths)
    if nrand_init>state_space_size
        nrand_init = state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
    state_ranges = [shuffle(sigma.state_range) for sigma in netstate]
    state_space = product(state_ranges...)
    i = 0
    for Xa in cycle(state_space)
        if rand()>1/length(netstate[end].state_range)
            continue
        end
        for (k,sigma) in enumerate(netstate)
            init_node!(sigma,Xa[k])
        end
        Xa_tag = join(map(string,round(Int,[Xa...])))
        if haskey(pairs[1],Xa_tag)
            continue
        end
        for t in 1:tau
            Xatp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
            for k in eachindex(netstate)
                init_node!(netstate[k],Xatp1[k])
            end
        end
        for H in 1:Hmax
            Xb = deepcopy([Xa...])
            for j in shuffle(index_nodes_to_flip)[1:H]
                Xb[j] = rand(filter(x -> x!=Xa[j],netstate[j].state_range))
            end
            Xb_tag = join(map(string,round(Int,Xb)))
            if in((Xa_tag => Xb_tag),pairs[H]) || in((Xb_tag => Xa_tag),pairs[H])
                nredundant += 1
            end
            pairs[H][Xa_tag] = Xb_tag
            for k in eachindex(netstate)
                init_node!(netstate[k],Xb[k])
            end
            for t in 1:tau
                Xbtp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
                for k in eachindex(netstate)
                    init_node!(netstate[k],Xbtp1[k])
                end
            end
            push!(Hmap,(H/N,sum(Xatp1.!=Xbtp1)/N))
            if mod(i,10000)==0
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
            end
        end
        i += 1
        if i == nrand_init
            break
        end
    end
    println("$nredundant/$(Hmax*nrand_init) configuration pairs redundant")
    return reshape(Hmap,Hmax,nrand_init)
end

function derrida_map3!(netstate::Vector{BNode},nrand_init::Integer=prod([length(sigma.state_range) for sigma in netstate]);tau::Integer=1,Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netstate]))
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,netstate)
    N = length(index_nodes_to_flip)
    1<=Hmax<=N || error("Hmax must be less or equal than network size")
    Hmap = Tuple{Float64,Float64}[]
    pairs = [Dict{AbstractString,AbstractString}() for i in 1:Hmax]
    # tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    nredundant = 0
    Xa = Vector{Integer}(length(netstate))
    Xb = Vector{Integer}(length(netstate))
    Xatp1 = Vector{Integer}(length(netstate))
    Xbtp1 = Vector{Integer}(length(netstate))
    range_lengths = [length(sigma.state_range) for sigma in netstate]
    state_space_size = prod(range_lengths)
    if nrand_init>state_space_size
        nrand_init = state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
    pat_div = reverse(round(Int,cumprod(reverse(range_lengths))./reverse(range_lengths)))
    state_space_number = shuffle(collect(1:state_space_size))
    i = 0
    for conf_i in state_space_number
        conf_indices = mod(div((conf_i-1),pat_div),range_lengths)
        Xa = [sigma.state_range[conf_indices[k]+1] for (k,sigma) in enumerate(netstate)]
        Xa_tag = join(map(string,round(Int,Xa)))
        for k in eachindex(netstate)
            init_node!(netstate[k],Xa[k])
        end
        for t in 1:tau
            Xatp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
            for k in eachindex(netstate)
                init_node!(netstate[k],Xatp1[k])
            end
        end
        for H in 1:Hmax
            Xb = deepcopy(Xa)
            for j in shuffle(index_nodes_to_flip)[1:H]
                Xb[j] = rand(filter(x -> x!=Xa[j],netstate[j].state_range))
            end
            Xb_tag = join(map(string,round(Int,Xb)))
            if in((Xa_tag => Xb_tag),pairs[H]) || in((Xb_tag => Xa_tag),pairs[H])
                nredundant += 1
            end
            pairs[H][Xa_tag] = Xb_tag
            for k in eachindex(netstate)
                init_node!(netstate[k],Xb[k])
            end
            for t in 1:tau
                Xbtp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
                for k in eachindex(netstate)
                    init_node!(netstate[k],Xbtp1[k])
                end
            end
            push!(Hmap,(H/N,sum(Xatp1.!=Xbtp1)/N))
            if mod(i,10000)==0
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
            end
        end
        i += 1
        if i == nrand_init
            break
        end
    end
    println("$nredundant/$(Hmax*nrand_init) configuration pairs redundant")
    return reshape(Hmap,Hmax,nrand_init)
end

function derrida_map_att3!(netstate::Vector{BNode},state_space_tags::Set{Int},nrand_init::Integer=length(state_space_tags);tau::Integer=1,Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netstate]))
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,netstate)
    N = length(index_nodes_to_flip)
    net_size = length(netstate)
    1<=Hmax<=N || error("Hmax must be less or equal than network size")
    Hmap = Tuple{Float64,Float64}[]
    pairs = [Dict{AbstractString,AbstractString}() for i in 1:Hmax]
    tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    nredundant = 0
    Xa = Vector{Integer}(length(netstate))
    Xb = Vector{Integer}(length(netstate))
    Xatp1 = Vector{Integer}(length(netstate))
    Xbtp1 = Vector{Integer}(length(netstate))
    sampled_confs = shuffle(collect(state_space_tags))
    if nrand_init>length(state_space_tags)
        nrand_init = length(state_space_tags)
    end
    i = 0
    for conf_i in sampled_confs
        Xa = reverse(digits(conf_i,tag_base,net_size))
        Xa_tag = join(map(string,round(Int,Xa)))
        for k in eachindex(netstate)
            init_node!(netstate[k],Xa[k])
        end
        Xatp1 = evolve_net!(netstate,tau,keep=false)
        for H in 1:Hmax
            belongs_to_att = false
            already_visitted = true
            while !belongs_to_att && already_visitted
                Xb = deepcopy(Xa)
                for j in sample(index_nodes_to_flip, H, replace = false, ordered = true)
                    Xb[j] = rand(filter(x -> x!=Xa[j],netstate[j].state_range))
                end
                Xb_tag = parse(Int,join(map(string,Xa)),tag_base)
                if in(Xb_tag,state_space_tags)
                    belongs_to_att = true
                end
                if in((Xa_tag => Xb_tag),pairs[H]) || in((Xb_tag => Xa_tag),pairs[H])
                    nredundant += 1
                    belongs_to_att = false
                else
                    pairs[H][Xa_tag] = Xb_tag
                    already_visitted = false
                    for k in eachindex(netstate)
                        init_node!(netstate[k],Xb[k])
                    end
                    Xbtp1 = evolve_net!(netstate,tau,keep=false)
                    push!(Hmap,(H/N,sum(Xatp1.!=Xbtp1)/N))
                    push!(Hmap,(H/N,sum(Xatp1.!=Xbtp1)/N))
                    if mod(i,10000)==0
                        @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
                    end
                end
            end
        end
        i += 1
        if i == nrand_init
            break
        end
    end
    println("$nredundant revisitted configuration pairs were discarded")
    return reshape(Hmap,Hmax.nrand_init)
end

function full_derrida_map!(netstate::Vector{BNode},nrand_init::Integer;tau::Integer=1,Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netstate]))
    nconf = prod([length(sigma.state_range) for sigma in netstate])
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,netstate)
    N = length(index_nodes_to_flip)
    1<=Hmax<=N || error("Hmax must be less or equal than network size")
    Hmap = Dict(i=>Tuple{AbstractFloat,AbstractFloat}[] for i in 1:Hmax)
    nredundant = 0
    Xatp1 = Vector{Integer}(length(netstate))
    Xbtp1 = Vector{Integer}(length(netstate))
    state_ranges = [shuffle(sigma.state_range) for sigma in netstate]
    state_space = collect(product(state_ranges...))
    i = 0
    for Xa in product([shuffle(sigma.state_range) for sigma in netstate]...)
        for (k,sigma) in enumerate(netstate)
            init_node!(sigma,Xa[k])
        end
        for t in 1:tau
            Xatp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
            for k in eachindex(netstate)
                init_node!(netstate[k],Xatp1[k])
            end
        end
        for Xb in product([shuffle(sigma.state_range) for sigma in netstate]...)
            if Xa==Xb
                continue
            end
            for (k,sigma) in enumerate(netstate)
                init_node!(sigma,Xb[k])
            end
            for t in 1:tau
                Xbtp1 = [evolve_node(netstate,nsigma) for nsigma in 1:length(netstate)]
                for k in eachindex(netstate)
                    init_node!(netstate[k],Xbtp1[k])
                end
            end
            H = sum(Xa.!=Xb)
            push!(Hmap[H],(H/N,sum(Xatp1.!=Xbtp1)/N))
            i+=1
            if mod(i,10000)==0
                @printf("%0.2f%% progress. H(0)=%d,%d\n", 100*i/nrand_init,H,length(Hmap[H]))
            end
            if i == nrand_init
                break
            end
        end
    end
    for H in sort(collect(keys(Hmap)))
        @show H,length(Hmap[H])
    end
    return Hmap
end


function derrida_map_mod!(netWTstate::Vector{BNode},netMutantstate::Vector{BNode},nrand_init::Integer=prod([length(sigma.state_range) for sigma in netMutantstate]);tau::Integer=1,fix_targets::Vector{Symbol}=Symbol[],Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netWTstate]))
    index_nodes_to_fix = Integer[]
    for id in fix_targets
        i = findfirst(x->x.id==string(id),netWTstate)
        j = findfirst(x->x.id==string(id),netMutantstate)
        (i != 0 && j != 0) || error("Trying to set unexisting $id node")
        i == j || error("Comparative derrida map requires networks to have the same size and same node order")
        push!(index_nodes_to_fix,i)
	end
    Hmin = length(index_nodes_to_fix)
    index_nodes_to_flip = find(x->length(x.state_range)>1, netWTstate)
    N = length(netWTstate)
    Hmin<=Hmax<=N || error("Hmax must be greater or equal than size of fixed targets, and equal or less than network size")
    deleteat!(index_nodes_to_flip,findin(index_nodes_to_flip,index_nodes_to_fix))
    Hmap = Tuple{Float64,Float64}[]
    size_Mutant_space_state = prod([length(sigma.state_range) for sigma in netMutantstate])
    if nrand_init>size_Mutant_space_state
        nrand_init = size_Mutant_space_state
    end
    pairs = Dict(i=>Dict{AbstractString,AbstractString}() for i in Hmin:Hmax)
    nredundant = 0
    XWTtp1 = Vector{Integer}(length(netWTstate))
    XMTtp1 = Vector{Integer}(length(netMutantstate))
    XWT = Vector{Integer}(length(netWTstate))
    XMT = Vector{Integer}(length(netMutantstate))
    i = 0
    while i < nrand_init
        XMT = rand_init_net!(netMutantstate)
        XMT_tag = join(map(string,round(Int,XMT)))
        for j in eachindex(netWTstate)
            init_node!(netWTstate[j],netMutantstate[j].state)
        end # the KO configuration is copied into the wild type (unperturbed) network
        for j in index_nodes_to_fix
            init_node!(netWTstate[j],rand(filter(x -> x!=netMutantstate[j].state,netWTstate[j].state_range)))
        end # those nodes to be shut off in the knock out network are set in the wild type network to any other state different than that
        tmpXWT = [sigma.state for sigma in netWTstate]
        for H in Hmin:Hmax
            already_visitted = true
            while already_visitted
                for j in eachindex(tmpXWT)
                    init_node!(netWTstate[j],tmpXWT[j])
                end # the KO configuration is copied into the wild type (unperturbed) network
                for j in sample(index_nodes_to_flip, H - Hmin, replace = false, ordered = true)
                    init_node!(netWTstate[j],rand(filter(x -> x!=netMutantstate[j].state,netWTstate[j].state_range)))
                end # the rest of the wildtype nodes are randomly set, up to H changes are done
                XWT = [sigma.state for sigma in netWTstate]
                XWT_tag = join(map(string,round(Int,XWT)))
                if (XMT_tag != XMT_tag) & (in((XWT_tag => XMT_tag),pairs[H]) || in((XMT_tag => XWT_tag),pairs[H]))
                    nredundant += 1
                else
                    pairs[H][XMT_tag] = XWT_tag
                    already_visitted = false
                end
            end
            for t in 1:tau
                XWTtp1 = [evolve_node(netWTstate,nsigma) for nsigma in 1:length(netWTstate)]
                for k in eachindex(netWTstate)
                    init_node!(netWTstate[k],XWTtp1[k])
                end
            end
            for t in 1:tau
                XMTtp1 = [evolve_node(netMutantstate,nsigma) for nsigma in 1:length(netMutantstate)]
                for k in eachindex(netMutantstate)
                    init_node!(netMutantstate[k],XMTtp1[k])
                end
            end
            # push!(Hmap,(H/N,sum(XWTtp1.!=XMTtp1)/N))
            push!(Hmap,(sum(XWT.!=XMT)/N,sum(XWTtp1.!=XMTtp1)/N))
            if mod(i,10000)==0
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
            end
        end
        i += 1
    end
    println("$nredundant revisitted configuration pairs were discarded")
    return reshape(Hmap,Hmax-Hmin+1,nrand_init)
end

function derrida_map_mod3!(netWTstate::Vector{BNode},netMutantstate::Vector{BNode},nrand_init::Integer=prod([length(sigma.state_range) for sigma in netMutantstate]);tau::Integer=1,fix_targets::Vector{Symbol}=Symbol[],Hmax::Integer=sum([length(sigma.state_range)>1 for sigma in netWTstate]))
    nconfMT = prod([length(sigma.state_range) for sigma in netMutantstate])
    if nrand_init>nconfMT
        nrand_init = nconfMT
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
    index_nodes_to_fix = Integer[]
    for id in fix_targets
        i = findfirst(x->x.id==string(id),netWTstate)
        j = findfirst(x->x.id==string(id),netMutantstate)
        (i != 0 && j != 0) || error("Trying to set unexisting $id node")
        i == j || error("Comparative derrida map requires networks to have the same size and same node order")
        push!(index_nodes_to_fix,i)
    end
    Hmin = length(index_nodes_to_fix)
    index_nodes_to_flip = find(x->length(x.state_range)>1, netWTstate)
    N = length(netWTstate)
    Hmin<=Hmax<=N || error("Hmax must be greater or equal than size of fixed targets, and equal or less than network size")
    deleteat!(index_nodes_to_flip,findin(index_nodes_to_flip,index_nodes_to_fix))
    Hmap = Tuple{Float64,Float64}[]
    pairs = Dict(i=>Dict{AbstractString,AbstractString}() for i in Hmin:Hmax)
    nredundant = 0
    XWTtp1 = Vector{Integer}(length(netWTstate))
    XMTtp1 = Vector{Integer}(length(netMutantstate))
    XWT = Vector{Integer}(length(netWTstate))
    XMT = Vector{Integer}(length(netMutantstate))
    range_lengths = [length(sigma.state_range) for sigma in netMutantstate]
    pat_div = reverse(round(Int,cumprod(reverse(range_lengths))./reverse(range_lengths)))
    state_space_number = shuffle(collect(1:prod(range_lengths)))
    i = 0
    for conf_i in state_space_number
        conf_indices = mod(div((conf_i-1),pat_div),range_lengths)
        XMT = [sigma.state_range[conf_indices[k]+1] for (k,sigma) in enumerate(netMutantstate)]
        XMT_tag = join(map(string,round(Int,XMT)))
        for j in eachindex(netMutantstate)
            init_node!(netWTstate[j],XMT[j])
        end # the KO configuration is copied into the wild type (unperturbed) network
        for j in index_nodes_to_fix
            init_node!(netWTstate[j],rand(filter(x -> x!=XMT[j],netWTstate[j].state_range)))
        end # those nodes to be shut off in the knock out network are set in the wild  type network to any other state different than that
        tmpXWT = [sigma.state for sigma in netWTstate]
        for H in Hmin:Hmax
            already_visitted = true
            for j in eachindex(netMutantstate)
                init_node!(netMutantstate[j],XMT[j])
            end # the KO configuration is copied into the wild type (unperturbed) network
            while already_visitted
                for j in eachindex(tmpXWT)
                    init_node!(netWTstate[j],tmpXWT[j])
                end # the KO configuration is copied into the wild type (unperturbed) network
                for j in sample(index_nodes_to_flip, H - Hmin, replace = false, ordered = true)
                    init_node!(netWTstate[j],rand(filter(x -> x!=XMT[j],netWTstate[j].state_range)))
                end # the rest of the wildtype nodes are randomly set, up to H changes are done
                XWT = [sigma.state for sigma in netWTstate]
                XWT_tag = join(map(string,round(Int,XWT)))
                if (XMT_tag != XMT_tag) & (in((XWT_tag => XMT_tag),pairs[H]) || in((XMT_tag => XWT_tag),pairs[H]))
                    nredundant += 1
                else
                    pairs[H][XMT_tag] = XWT_tag
                    already_visitted = false
                end
            end
            for t in 1:tau
                XWTtp1 = [evolve_node(netWTstate,nsigma) for nsigma in 1:length(netWTstate)]
                for k in eachindex(netWTstate)
                    init_node!(netWTstate[k],XWTtp1[k])
                end
            end
            for t in 1:tau
                XMTtp1 = [evolve_node(netMutantstate,nsigma) for nsigma in 1:length(netMutantstate)]
                for k in eachindex(netMutantstate)
                    init_node!(netMutantstate[k],XMTtp1[k])
                end
            end
            # push!(Hmap,(H/N,sum(XWTtp1.!=XMTtp1)/N))
            push!(Hmap,(sum(XWT.!=XMT)/N,sum(XWTtp1.!=XMTtp1)/N))
            if mod(i,10000)==0
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
            end
        end
        i += 1
        if i==nrand_init
            break
        end
    end
    println("$nredundant revisitted configuration pairs were discarded")
    return reshape(Hmap,Hmax-Hmin+1,nrand_init)
end


# get_basins explores all possible configurations and returns a transition graph, with a dictionary storing labels for
# each configuration

function get_basins!(netstate::Vector{BNode},out_dir::AbstractString,root_name::AbstractString;trav::Symbol=:rw,update_mode::Symbol=:sy,max_traj::Integer=maximum(vcat([sigma.tau_in for sigma in netstate]...)))
    trav == :rw || trav == :l || error("Invalid option for traversal method")
    net_size = length(netstate)
    state_ranges = [sigma.state_range for sigma in netstate]
    itstspace = product(state_ranges...)
    size_space_state = length(itstspace)
    tags_nid = update_mode == :sy ? Dict{Integer,Integer}() : Dict{Vector{Integer},Integer}()
    nid_tags = update_mode == :sy ? Dict{Integer,Integer}() : Dict{Integer,Vector{Integer}}()
    tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    transition_graph = DiGraph()
    println("Calculating attraction landscape for $root_name")
    tau_max = maximum(vcat([sigma.tau_in for sigma in netstate]...))
    net_hist_size = tau_max + 1
    if update_mode == :sy
        for leave in itstspace
            for k in eachindex(netstate)
                init_node!(netstate[k],leave[k])
            end
            leave_tag = parse(Int,string(map(x->round(Int,x),leave)...),tag_base)
            shoot_tag = parse(Int,string(round(Int,shoot)...),tag_base)
            if !haskey(tags_nid,leave_tag)
                add_vertex!(transition_graph)
                nvertex = nv(transition_graph)
                nid_tags[nvertex] = leave_tag
                tags_nid[leave_tag] = nvertex
            end
            shoot = evolve_net!(netstate,1,keep=false,update_mode=update_mode)
            if !haskey(tags_nid,shoot_tag)
                add_vertex!(transition_graph)
                nvertex = nv(transition_graph)
                nid_tags[nvertex] = shoot_tag
                tags_nid[shoot_tag] = nvertex
            end
            add_edge!(transition_graph,tags_nid[leave_tag],tags_nid[shoot_tag])
            if mod(nv(transition_graph),100000)==0
                @printf "%0.2f%% configurations explored\n" 100*nv(transition_graph)/size_space_state
            end
        end
    else
        hX = Array{Int8}(length(netstate),net_hist_size)
        hX_tag = Vector{Int}(net_hist_size)
        hXtp1 = Array{Int8}(length(netstate),net_hist_size)
        hXtp1_tag = Vector{Int}(net_hist_size)
        time_steps = 1
        sim_buffer = Array{Int8}(length(netstate),time_steps + net_hist_size)
        for (i,X) in enumerate(itstspace)
            hX[:] = init_net!(netstate;init_t0=:default,map((x,y)->Symbol(x)=>y,(sigma.id for sigma in netstate),X)...)
            # for (nsigma,sigma) in enumerate(netstate)
            #     init_node!(sigma,X[nsigma])
            #     hX[nsigma,:] = sigma.state
            # end
            for t in 1:net_hist_size
                hX_tag[t] = parse(Int,string(round(Int,hX[:,t])...),tag_base)
            end
            if !haskey(tags_nid,hX_tag)
                add_vertex!(transition_graph)
                nvertex = nv(transition_graph)
                nid_tags[nvertex] = hX_tag
                tags_nid[hX_tag] = nvertex
            end
            # hXtp1_tag[:] = hX_tag[:]
            # shift!(hXtp1_tag)
            for t in 0:max_traj-1
                hXtp1[:] = evolve_net!(netstate,1,keep=false,update_mode=update_mode,sim=sim_buffer,t0=t)
                for ti in 1:net_hist_size
                    hXtp1_tag[ti] = parse(Int,string(round(Int,hXtp1[:,ti])...),tag_base)
                end
                # Xtp1_tag = parse(Int,string(round(Int,Xtp1)...),tag_base)
                # push!(hXtp1_tag,Xtp1_tag)
                if !haskey(tags_nid,hXtp1_tag)
                    add_vertex!(transition_graph)
                    nvertex = nv(transition_graph)
                    nid_tags[nvertex] = hXtp1_tag
                    tags_nid[hXtp1_tag] = nvertex
                end
                add_edge!(transition_graph,tags_nid[hX_tag],tags_nid[hXtp1_tag]) || break
                hX_tag[:] = hXtp1_tag[:]
                # shift!(hXtp1_tag)
            end
            if mod(i,100000)==0
                @printf "%0.2f%% initial conditions explored\n" 100*i/size_space_state
            end
        end
    end
    size_space_state = nv(transition_graph)
    println("Decomposing transition graph into attraction basins")
    basin_node_list = weakly_connected_components(transition_graph)
    Natt = length(basin_node_list)
    println(Natt," attractors found")
    sort!(basin_node_list,by=length,rev=true)
    println("Saving individual basins")
    if !isdir(out_dir)
        mkpath(out_dir)
    end
    g = GZip.open(joinpath(out_dir,"$(root_name)_att_tagging.dat.gz"),"w")
    h = GZip.open(joinpath(out_dir,"$(root_name)_trans_dist.dat.gz"),"w")
    try
        basin_tree_list = Vector{LightGraphs.DiGraph}(Natt)
        basin_tags = Vector{Vector{Integer}}(Natt)
        att_list = Vector{Vector{Integer}}(Natt)
        periods = Vector{Integer}(Natt)
        branch_sizes = Vector{Vector{Integer}}(Natt)
        patt_basins = Vector{Vector{Integer}}(Natt)
        branch_tree_list = Vector{Vector{LightGraphs.DiGraph}}(Natt)
        branch_patt = Vector{Vector{Integer}}(Natt)
        for (bi,basin) in enumerate(basin_node_list)
            basin_tree_list[bi] = induced_subgraph(transition_graph,basin)[1]
            save(joinpath(out_dir,"$(root_name)_att_$(bi).gml.gz"),basin_tree_list[bi],join([root_name, bi],"_"),:gml,compress = true)
            println("$(root_name)_att_$(bi) basin saved")
            if update_mode == :sy
                writedlm(g,hcat(fill(bi,nv(basin_tree_list[bi])),[nid_tags[v] for v in basin]))
            else
                writedlm(g,hcat(fill(bi,nv(basin_tree_list[bi])),reshape([nid_tags[v][d] for d in 1:net_hist_size for v in basin],length(basin),net_hist_size)))
            end
            println("Computing attractor of basin $bi")
            probe = basin[1]
            path_to_att = saw(transition_graph,probe,length(basin))
            path_in_att = saw(transition_graph,path_to_att[end],length(basin))
            att_list[bi] = path_in_att
            periods[bi] = length(path_in_att)
            println("Attractor $bi has period ",periods[bi])
            patt_basin = Vector{Integer}()
            for point_in_att in path_in_att
                push!(patt_basin, findfirst(basin,point_in_att))
            end
            patt_basins[bi] = patt_basin
        end
        transition_graph = DiGraph()
        for (bi,basin_tree) in enumerate(basin_tree_list)
            println("Calculating branches associated to attractor $bi")
            edges_att = hcat(circshift(patt_basins[bi],1),patt_basins[bi])
            for pi in 1:length(patt_basins[bi])
                rem_edge!(basin_tree,edges_att[pi,:]...)
            end
            branch_node_list = weakly_connected_components(basin_tree)
            branch_tree_list[bi] = Vector{LightGraphs.DiGraph}(periods[bi])
            patt_to_look = deepcopy(patt_basins[bi])
            branch_patt[bi] = Vector{Integer}(periods[bi])
            for (bri,branch) in enumerate(branch_node_list)
                for point_in_att in patt_to_look
                    if point_in_att in branch
                        branch_patt[bi][bri] = findfirst(patt_basins[bi],point_in_att)
                        branch_tree_list[bi][bri] = induced_subgraph(basin_tree,branch)[1]
                        # println("Branch $bi:$bri has $(nv(branch_tree_list[bi][bri])) nodes")
                        deleteat!(patt_to_look,findfirst(patt_to_look,point_in_att))
                        break
                    end
                end
            end
        end
        basin_tree_list = Vector{LightGraphs.DiGraph}()
        for (bi,basin_tree) in enumerate(branch_tree_list)
            println("Computing transients distribution associated to attractor $bi")
            dists = Integer[]
            brsize = Vector{Integer}(periods[bi])
            for (bri,branch_tree) in enumerate(basin_tree)
                println("Measuring transient distances for branch $bi:$bri")
                dists_br = Integer[]
                if (trav==:rw)
                    maxbr_size = nv(branch_tree)
                    for pbi in vertices(branch_tree)
                      # push!(dists_br,length(a_star(basin_tree,pbranch,branch_patt[bri])))
                        push!(dists_br,length(non_backtracking_randomwalk(branch_tree,pbi,maxbr_size)-1))
                        if mod(pbi,100000)==0
                            @printf("%0.2f%% of branch %d:%d explored\n",100*pbi/maxbr_size, bri, bi)
                        end
                    end
                    filter!(x->x!=0,dists_br)
                else
                    m = mapreduce(x->in_neighbors(branch_tree,x),vcat,vertices(branch_tree))
                    d = sparsevec(m,ones(Integer,length(m)),nv(branch_tree))
                    l = 1
                    println("Level $l")
                    while !isempty(m)
                        m = mapreduce(x->in_neighbors(branch_tree,x),vcat,m)
                        d[m] += 1
                        l+=1
                        println("Level $l")
                    end
                    dists_br = nonzeros(d)
                end
                append!(dists,dists_br)
                brsize[branch_patt[bi][bri]] = length(dists_br)
            end
            branch_sizes[bi] = brsize
            writedlm(h,dists',',')
            println("Maximum transient distance:$(maximum(dists))")
        end
    finally
        close(g)
        close(h)
    end
    f = open(joinpath(out_dir,"$(root_name)_att_stats.txt"),"w")
    sum_att = sum(periods)
    sum_basin = size_space_state - sum_att
    try
        println("Writing statistics")
        for bi in eachindex(att_list)
            if bi > 1
                write(f,"\n")
            end
            write(f,"Attractor:$bi\t")
            period = periods[bi]
            write(f,"Period:$period\t")
            basin_size = mapreduce(x->nv(x),+,branch_tree_list[bi])-period
            write(f,"Basin size:$basin_size/$sum_basin ($(format(basin_size/sum_basin*100.0,precision=2,zeropadding=true))%)\n")
            if update_mode==:sy
                for (pi,conf_id) in enumerate(att_list[bi])
                    conf = reverse(digits(nid_tags[conf_id],tag_base,net_size))
                    write(f,"$(join(conf,","))\t$(branch_sizes[bi][pi])\n")
                end
            else
                for (pi,conf_id) in enumerate(att_list[bi])
                    conf = reverse(digits(nid_tags[conf_id][end],tag_base,net_size))
                    write(f,"$(join(conf,","))\t$(branch_sizes[bi][pi])\n")
                end
            end
        end
    finally
        close(f)
    end
    # return basin_tags,basin_tree_list,att_list
    return transition_graph
end

# get_attractor finds attractors exploring multiple random initial conditions (nrand_init) within a simulation
# time range limited by max_steps, after a transient given by transient_size
function find_attractor!(netstate::Vector{BNode},transient_size::Integer,max_steps::Integer,nrand_init::Integer;update_mode::Symbol=:sy)
    i = 0
    tag_base = maximum([length(sigma.state_range) for sigma in netstate])
    tau_max = update_mode==:sy ? 0 : maximum(vcat([sigma.tau_in for sigma in netstate]...))
    net_hist_size = tau_max + 1
    net_size = length(netstate)
    att_list = Vector{Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}}()
    hX = Array{Int8}(net_size,net_hist_size)
    sim_trans_buffer = Array{Int8}(net_size,transient_size + net_hist_size)
    time_steps = 1
    sim_att_buffer = Array{Int8}(net_size,time_steps + net_hist_size)
    unsuccess = 0
    while i < nrand_init
        init_net!(netstate;init_t0=:rand,init_history=:rand)
        p_att = Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}()
        hX[:] = evolve_net!(netstate,transient_size,keep=false,update_mode=update_mode,sim=sim_trans_buffer)
        probe = tuple(map(x->tuple(x...),(hX[:,ti] for ti in 1:net_hist_size))...)
        push!(p_att,probe)
        add_attractor = false
        already_found = false
        attractor_found = false
        tmax = 1
        for t in transient_size:transient_size+max_steps-1
            hX[:] = evolve_net!(netstate,1,keep=false,update_mode=update_mode,t0=t,sim=sim_att_buffer)
            ntp1 = tuple(map(x->tuple(x...),(hX[:,ti] for ti in 1:net_hist_size))...)
            tmax = t
            if probe[end] == ntp1[end]
                # println("Likely attractor found")
                for att in att_list
                    for p_att_i in p_att
                        if p_att_i in att
                            already_found = true
                            break
                        end
                    end
                    if already_found
                        break
                    end
                end
                attractor_found = true
                if !already_found
                    push!(att_list,p_att)
                    println("Attractor $(size(att_list,1)) has period $(size(p_att,1))")
                    for (j,att_i) in enumerate(p_att)
                        println(join(att_i[end],","))
                    end
                end
                break
            end
            push!(p_att,ntp1)
        end
        i += 1
        if !(already_found||attractor_found)
            println(tmax)
            print("probe:")
            println(join(probe[end],","))
            print("ntp1 :")
            println(join(hX[:,end],","))
            #  println("No attractors after reaching maxtimum time steps")
            unsuccess += 1
            # println("No attractor found for initial condition $i:")
            # for probe_i in probe
            #     conf = reverse(digits(probe_i,tag_base,net_size))
            #     println(join(conf,","))
            # end
        end
    end
    println(unsuccess," unsucessful tries after reaching maximmum time steps set")
    return att_list
end
