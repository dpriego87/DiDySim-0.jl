import GZip
using StatsBase
using IterTools
using Combinatorics
using LightGraphs
using Formatting
#derrida_map calculates derrida map for a net up to a  H hamming distance
function derrida_map!{T<:state_type}(net::Net2{T},nrand_init::Int=prod([length(sigma.state_range) for sigma in net.nodes]);tau::Int=1,update_mode::Symbol=:sy,init_h::Symbol=:rand,Hmax::Int=sum(length(sigma.state_range)>1 for sigma in net.nodes)*(1 + (update_mode==:sy ? 0 : maximum(net.in_taus))))
    index_nodes_to_flip = find(sigma->length(sigma.state_range)>1,net.nodes)
    net_hist_size = 1 + (update_mode==:sy ? 0 : maximum(net.in_taus))
    Nmax = length(index_nodes_to_flip)*net_hist_size
    net_size = length(net.nodes)
    1<=Hmax<=Nmax || error("Hmax must be less or equal than network size")
    state_space_size = prod([length(sigma.state_range) for sigma in net.nodes])
    if nrand_init>state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
        nrand_init = state_space_size
    end
    mut_indices = collect(product(index_nodes_to_flip,1:net_hist_size))
    hXa = Array{Int}(net_size,net_hist_size)
    hXatp = Array{Int}(net_size,net_hist_size)
    hXb = Array{Int}(net_size,net_hist_size)
    hXbtp = Array{Int}(net_size,net_hist_size)
    # sim_buffer = Array{T}(net_size,tau + net_hist_size)
    Hmap = Tuple{Float64,Float64}[]
    i = 0
    while i < nrand_init
        rand_init_net!(net;update_mode=update_mode,init_h=init_h)
        t = rand(1:net_hist_size)-1
        evolve_net!(net,net_hist_size-1,keep=false,update_mode=update_mode,t_init=t)
        copy!(hXa,net.sim_buffer[:,end-net_hist_size+1:end])
        evolve_net!(net,tau,keep=false,update_mode=update_mode,t_init=t)
        copy!(hXatp,net.sim_buffer[:,end-net_hist_size+1:end])
        for H in 1:Hmax
            copy!(hXb,hXa)
            for (n,l) in sample(mut_indices, H, replace = false)
                hXb[n,l] = rand(filter(x -> x!=hXa[n,l],net.nodes[n].state_range))
            end
            # print("H=$H :",sum(hXa.!=hXb)," -> ")
            # init_net!(net;update_mode=update_mode,init_array=vec(hXb))
            for ni in 1:net_size
                net.sim_buffer[ni,end-net_hist_size+1:end] = collect(hXb[ni,:])
            end
            evolve_net!(net,tau,keep=false,update_mode=update_mode,t_init=t)
            copy!(hXbtp,net.sim_buffer[:,end-net_hist_size+1:end])
            # println(sum(net.sim_buffer[:,end-net_hist_size+1:end].!=hXatp))
            push!(Hmap,(H/Nmax,sum(hXatp.!=hXbtp)/Nmax))
            if mod(i,10000)==1
                @printf("%0.2f%% progress, %d out of %d random initial conditions explored at distance H(0)=%d\n", 100*i/nrand_init, i, nrand_init,H)
                println("H(0)=$H:")
                # println(join(map(string,hXa),","))
                # println(join(map(string,hXb),","))
                println("H($tau)=",sum(hXatp.!=hXbtp),":")
                # println(join(map(string,hXatp),","))
                # println(join(map(string,hXbtp),","))
            end
        end
        i += 1
    end
    return reshape(Hmap,Hmax,nrand_init)
end

function derrida_map_mod!{T<:state_type}(netWT::Net2{T},nrand_init::Int=prod([length(sigma.state_range) for sigma in netWT.nodes]);mutants::Vector{Pair{String,Symbol}} = Vector{Pair{String,Symbol}}(),tau::Int=1,fix_targets::Bool=true,update_mode::Symbol=:sy,init_h::Symbol=:rand,Hmax::Int=sum([length(sigma.state_range)>1 for sigma in netWT.nodes]))
    netMT = deepcopy(netWT)
    for (id,pert) in mutants
        if pert == :KO
            knock_out!(netMT,id)
        elseif pert == :KD
            knock_down!(netMT,id)
        else
            error("Invalid method")
        end

    end
    Mutant_state_space_size = prod([length(sigma.state_range) for sigma in netMT.nodes])
    Mutant_state_space_size = length(combinations(1:Mutant_state_space_size,2))
    if nrand_init>Mutant_state_space_size
        nrand_init = Mutant_state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
    net_hist_size = 1 + (update_mode==:sy ? 0 : maximum(netWT.in_taus))
    index_nodes_to_fix = Integer[]
    if fix_targets
        for (id,_) in mutants
            i = findfirst(x->x.id==id,netWT.nodes)
            j = findfirst(x->x.id==id,netMT.nodes)
            (i != 0 && j != 0) || error("Trying to set unexisting $id node")
            #i == j || error("Comparative derrida map requires networks to have the same size and same node order")
            push!(index_nodes_to_fix,i)
        end
    end
    Hmin = length(index_nodes_to_fix)
    index_nodes_to_flip = find(x->length(x.state_range)>1, netWT.nodes)
    deleteat!(index_nodes_to_flip,findin(index_nodes_to_flip,index_nodes_to_fix))
    NWT = length(netWT.nodes)
    NMT = length(netMT.nodes)
    net_size = NWT
    #Nmax = length(index_nodes_to_flip)*net_hist_size
    Nmax = NWT*net_hist_size
    Hmin<=Hmax<=Nmax || error("Hmax must be greater or equal than size of fixed targets, and equal or less than network size")
    pairs = Dict(i=>Dict{AbstractString,AbstractString}() for i in Hmin:Hmax)
    mut_indices = shuffle(collect(product(index_nodes_to_flip,1:net_hist_size)))
    hXWT = Array{Int}(NWT,net_hist_size)
    hXWTtptau = collect(Array{Int}(NWT,net_hist_size) for _ in 1:tau)
    hXMT = Array{Int}(NMT,net_hist_size)
    hXMTtptau = collect(Array{Int}(NMT,net_hist_size) for _ in 1:tau)
    Hmap = Tuple{Int64,Vararg{Float64,2}}[]
    i = 0
    while i < nrand_init
        rand_init_net!(netMT;update_mode=update_mode,init_h=init_h)
        copy!(hXMT,netMT.sim_buffer[:,end-net_hist_size+1:end])
        t = rand(1:net_hist_size)-1
        netMTtimeseries = evolve_net!(netMT,tau,keep=true,update_mode=update_mode,t_init=t)
        for ti in 1:tau
            copy!(hXMTtptau[ti],netMTtimeseries[:,1+tau-net_hist_size+1:1+tau])
        end
        for H in Hmin:Hmax
            copy!(hXWT,hXMT)
            for (n,l) in sample(mut_indices, H-Hmin, replace = false)
                hXWT[n,l] = rand(filter(x -> x!=hXMT[n,l],netWT.nodes[n].state_range))
            end
            for n in index_nodes_to_fix
                hXWT[n,end] = rand(filter(x -> x!=hXMT[n,end],netWT.nodes[n].state_range))
            end # those nodes to be shut off in the knock out network are set in the wild type network to any other state different than that
            for ni in 1:NWT
                netWT.sim_buffer[ni,end-net_hist_size+1:end] = collect(hXWT[ni,:])
            end
            init_net!(netWT;t0_mode=:custom,update_mode=update_mode,init_h=:last)
            netWTtimeseries = evolve_net!(netWT,tau,keep=true,update_mode=update_mode,t_init=t)
            for ti in 1:tau
                copy!(hXWTtptau[ti],netWTtimeseries[:,1+tau-net_hist_size+1:1+tau])
                push!(Hmap,(ti,sum(hXWT.!=hXMT)/Nmax,sum(hXWTtptau[ti].!=hXMTtptau[ti])/Nmax))
            end
        end
        i += 1
        if mod(i,10000)==0
            @printf("%0.2f%% progress, %d out of %d random initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
        end
    end
    return reshape(Hmap,tau*(Hmax-Hmin+1),nrand_init)
end


# get_basins explores all possible configurations and returns a transition graph, with a dictionary storing labels for
# each configuration

function get_basins!(net::Net2,out_dir::AbstractString,root_name::AbstractString;trav::Symbol=:rw,update_mode::Symbol=:sy,max_traj::Integer=1 + (update_mode==:sy ? 0 : maximum(net.in_taus)))
    trav == :rw || trav == :l || error("Invalid option for traversal method")
    net_size = nv(net.graph)
    state_space_size = prod(length(sigma.state_range) for sigma in net.nodes)
    if state_space_size==0 && all(x->isinteger(x)&&x>0,length(sigma.state_range) for sigma in net.nodes)
        error("Full exploration of state space exceeds memory limit")
    end
    itstspace = product((sigma.state_range for sigma in net.nodes)...)
    state_space_size = length(itstspace)
    transition_graph = DiGraph()
    tags_nid = update_mode == :sy ? Dict{Int,Int}() : Dict{Vector{Int},Int}()
    nid_tags = update_mode == :sy ? Vector{Int}() : Vector{Vector{Int}}()
    tag_base = maximum(length(sigma.state_range) for sigma in net.nodes)
    println("Calculating attraction landscape for $root_name")
    tau_max = maximum(net.in_taus)
    net_hist_size = tau_max + 1
    if update_mode == :sy
        i = 0
        for leave in itstspace
            init_net!(net;update_mode=update_mode,init_array=round(Int,collect(leave)))
            leave_tag = parse(Int,string(map(x->round(Int,x),leave)...),tag_base)
            if !haskey(tags_nid,leave_tag)
                add_vertex!(transition_graph)
                i = nv(transition_graph)
                push!(nid_tags,leave_tag)
                tags_nid[leave_tag] = i
            end
            shoot = evolve_net!(net,1;keep=false,update_mode=update_mode)
            shoot_tag = parse(Int,string(round(Int,shoot)...),tag_base)
            if !haskey(tags_nid,shoot_tag)
                add_vertex!(transition_graph)
                i = nv(transition_graph)
                push!(nid_tags,shoot_tag)
                tags_nid[shoot_tag] = i
            end
            add_edge!(transition_graph,tags_nid[leave_tag],tags_nid[shoot_tag])
            if mod(i,100000)==0
                @printf "%0.3f%% configurations explored\n" 100*i/state_space_size
            end
        end
    elseif update_mode == :fd || update_mode == :md
        hX = Array{Int8}(net_size,net_hist_size)
        hX_tag = Vector{Int}(net_hist_size)
        hXtp1 = Array{Int8}(net_size,net_hist_size)
        hXtp1_tag = Vector{Int}(net_hist_size)
        time_steps = 1
        sim_buffer = Array{Int8}(net_size,time_steps + net_hist_size)
        for (i,X) in enumerate(itstspace)
            init_net!(net;init_h=:default,update_mode=update_mode,init_array=round(Int,collect(X)))
            hX[:] = net.sim_buffer
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
                push!(nid_tags, collect(hX_tag))
                tags_nid[collect(hX_tag)] = nvertex
                # println(hX_tag, " tag added, vertex ",tags_nid[hX_tag])
            end
            # hXtp1_tag[:] = hX_tag[:]
            # shift!(hXtp1_tag)
            for t in 1:max_traj
                hXtp1[:] = evolve_net!(net,1,keep=false,update_mode=update_mode,sim=sim_buffer,t_init=t-1)
                for ti in 1:net_hist_size
                    hXtp1_tag[ti] = parse(Int,string(round(Int,hXtp1[:,ti])...),tag_base)
                end
                # Xtp1_tag = parse(Int,string(round(Int,Xtp1)...),tag_base)
                # push!(hXtp1_tag,Xtp1_tag)
                if !haskey(tags_nid,hXtp1_tag)
                    add_vertex!(transition_graph)
                    nvertex = nv(transition_graph)
                    push!(nid_tags,collect(hXtp1_tag))
                    tags_nid[collect(hXtp1_tag)] = nvertex
                    # println(hXtp1_tag, " tag added, vertex ",tags_nid[hXtp1_tag])
                end
                # foreach(x->println(x," ",typeof(x)," ",hX_tag==x),keys(tags_nid))
                add_edge!(transition_graph,tags_nid[hX_tag],tags_nid[hXtp1_tag]) || break
                copy!(hX_tag,hXtp1_tag)
                # shift!(hXtp1_tag)
            end
            if mod(i,100000)==0
                @printf "%0.2f%% initial conditions explored\n" 100*i/state_space_size
            end
        end
    elseif update_mode == :st || update_mode == :mst
        hX = Array{Int8}(net_size,net_hist_size)
        hX_tag = Vector{Int}(net_hist_size+1)
        hXtp1 = Array{Int8}(net_size,net_hist_size)
        hXtp1_tag = Vector{Int}(net_hist_size+1)
        time_steps = 1
        sim_buffer = Array{Int8}(net_size,time_steps + net_hist_size)
        for (i,X) in enumerate(itstspace)
            init_net!(net;init_h=:default,update_mode=update_mode,init_array=round(Int,collect(X)))
            hX[:] = net.sim_buffer
            # for (nsigma,sigma) in enumerate(netstate)
            #     init_node!(sigma,X[nsigma])
            #     hX[nsigma,:] = sigma.state
            # end
            for t in 1:net_hist_size
                hX_tag[1+t] = parse(Int,string(round(Int,hX[:,t])...),tag_base)
            end
            hX_tag[1] = 0
            if !haskey(tags_nid,hX_tag)
                add_vertex!(transition_graph)
                nvertex = nv(transition_graph)
                push!(nid_tags, collect(hX_tag))
                tags_nid[collect(hX_tag)] = nvertex
                # println(hX_tag, " tag added, vertex ",tags_nid[hX_tag])
            end
            # hXtp1_tag[:] = hX_tag[:]
            # shift!(hXtp1_tag)
            # tmax = 0
            for t in 1:max_traj
                hXtp1[:] = evolve_net!(net,1,keep=false,update_mode=update_mode,sim=sim_buffer,t_init=t-1)
                for ti in 1:net_hist_size
                    hXtp1_tag[1+ti] = parse(Int,string(round(Int,hXtp1[:,ti])...),tag_base)
                end
                hXtp1_tag[1] = mod(t,net_hist_size)
                # Xtp1_tag = parse(Int,string(round(Int,Xtp1)...),tag_base)
                # push!(hXtp1_tag,Xtp1_tag)
                # println(t, " ",keys(tags_nid))
                if !haskey(tags_nid,hXtp1_tag)
                    add_vertex!(transition_graph)
                    nvertex = nv(transition_graph)
                    push!(nid_tags, collect(hXtp1_tag))
                    tags_nid[collect(hXtp1_tag)] = nvertex
                    # println(hXtp1_tag, " tag added, vertex ",tags_nid[hXtp1_tag])
                end
                # println(t, " ",keys(tags_nid))
                # tmax = t
                add_edge!(transition_graph,tags_nid[hX_tag],tags_nid[hXtp1_tag]) || break
                copy!(hX_tag,hXtp1_tag)
                # shift!(hXtp1_tag)
            end
            # println(tmax)
            if mod(i,100000)==0
                @printf "%0.2f%% initial conditions explored\n" 100*i/state_space_size
            end
        end
    end
    state_space_size = nv(transition_graph)
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
    basin_tree_list = Vector{LightGraphs.DiGraph}(Natt)
    basin_tags = Vector{Vector{Integer}}(Natt)
    att_list = Vector{Vector{Integer}}(Natt)
    periods = Vector{Integer}(Natt)
    branch_sizes = Vector{Vector{Integer}}(Natt)
    patt_basins = Vector{Vector{Integer}}(Natt)
    branch_tree_list = Vector{Vector{LightGraphs.DiGraph}}(Natt)
    branch_patt = Vector{Vector{Integer}}(Natt)
    try
        for (bi,basin) in enumerate(basin_node_list)
            basin_tree_list[bi] = induced_subgraph(transition_graph,basin)[1]
            save(joinpath(out_dir,"$(root_name)_att_$(bi).gml.gz"),basin_tree_list[bi],join([root_name, bi],"_"),:gml,compress = true)
            println("$(root_name)_att_$(bi) basin saved")
            if update_mode == :sy
                writedlm(g,hcat(fill(bi,nv(basin_tree_list[bi])),nid_tags[basin]))
            else
                writedlm(g,hcat(fill(bi,nv(basin_tree_list[bi])),reshape([nid_tags[v][d] for d in eachindex(nid_tags[1]) for v in basin],length(basin),length(nid_tags[1]))))
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
                            @printf("%0.2f%% of branch %d:%d explored\n",100*pbi/maxbr_size, bi, bri)
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
    sum_basin = state_space_size - sum_att
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
    return basin_tags,basin_tree_list,att_list
    # return transition_graph
end

# get_attractor finds attractors exploring multiple random initial conditions (nrand_init) within a simulation
# time range limited by max_steps, after a transient given by transient_size
function find_attractor!{T<:state_type}(net::Net2{T},transient_size::Integer,max_steps::Integer,nrand_init::Integer;update_mode::Symbol=:sy,forcing_rhythms::Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}=Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}(),init_constraints...)
    i = 0
    # tag_base = maximum([length(sigma.state_range) for sigma in net.nodes])
    net_hist_size = 1 + (update_mode==:sy ? 0 : maximum(net.in_taus))
    net_size = length(net.nodes)
    att_list = Vector{Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}}()
    hX = Array{T}(net_size,net_hist_size)
    sim_trans_buffer = Array{T}(net_size,transient_size + net_hist_size)
    time_steps = 1
    sim_att_buffer = Array{T}(net_size,time_steps + net_hist_size)
    unsuccess = 0
    basin_sizes = Int[]
    sigma_ext_forcing = collect(sigma_symbol=>rhythm_values[2] for (sigma_symbol,rhythm_values) in forcing_rhythms)
	init_constraints = vcat(init_constraints,sigma_ext_forcing)
    if update_mode in [:roa,:ga]
        error("Not yet implemented for random asyncrhonous updating")
    end
    while i < nrand_init
        rand_init_net!(net;update_mode=update_mode,init_h=:rand,init_constraints...)
        init_probe = tuple(map(x->tuple(x...),(net.sim_buffer[:,ti] for ti in 1:net_hist_size))...)
        p_att = Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}()
        hX[:] = evolve_net!(net,transient_size,keep=false,update_mode=update_mode,sim=sim_trans_buffer,forcing_rhythms=forcing_rhythms)
        probe = tuple(map(x->tuple(x...),(hX[:,ti] for ti in 1:net_hist_size))...)
        push!(p_att,probe)
        add_attractor = false
        already_found = false
        attractor_found = false
        tmax = 1
        for t in transient_size:transient_size+max_steps-1
            hX[:] = evolve_net!(net,1,keep=false,update_mode=update_mode,sim=sim_att_buffer,t_init=t,forcing_rhythms=forcing_rhythms)
            ntp1 = tuple(map(x->tuple(x...),(hX[:,ti] for ti in 1:net_hist_size))...)
            tmax = t
            if probe[end] == ntp1[end]
                # println("Likely attractor found")
                for (pi,att) in enumerate(att_list)
                    for p_att_i in p_att
                        if p_att_i in att
                            already_found = true
                            break
                        end
                    end
                    if already_found
                        if !(init_probe in p_att)
                            basin_sizes[pi]+=1
                        end
                        break
                    end
                end
                attractor_found = true
                if !already_found
                    push!(basin_sizes,1)
                    push!(att_list,p_att)
                    println("Attractor ",size(att_list,1)," has period ",size(p_att,1))
                    for att_i in p_att
                        println(join(att_i[end],","))
                    end
                end
                break
            end
            push!(p_att,ntp1)
        end
        i += 1
        if !(already_found||attractor_found)
            println("No attractor found After $tmax steps")
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
    if length(att_list)!=0
        println("Estimated basin sizes:")
        for (pi,att) in enumerate(att_list)
            println("Attractor $pi:",format(100.0*basin_sizes[pi]/sum(basin_sizes),precision=2,zeropadding=true),"%")
        end
    end
    println(unsuccess," unsuccessful tries after reaching maximmum time steps set")
    return att_list
end

function find_attractor2!{T<:state_type}(net::Net2{T},max_steps::Integer,nrand_init::Integer;update_mode::Symbol=:sy,forcing_rhythms::Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}=Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}(),init_constraints...)
    i = 0
    # tag_base = maximum([length(sigma.state_range) for sigma in net.nodes])
    net_hist_size = 1 + (update_mode==:sy ? 0 : maximum(net.in_taus))
    net_size = length(net.nodes)
    att_list = Vector{Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}}()
    hX = Array{T}(net_size,net_hist_size)
    sim_att_buffer = Array{T}(net_size,1 + net_hist_size)
    unsuccess = 0
    basin_sizes = Int[]
    transients = Vector{Vector{Int}}()
    sigma_ext_forcing = collect(sigma_symbol=>rhythm_values[2] for (sigma_symbol,rhythm_values) in forcing_rhythms)
	init_constraints = vcat(init_constraints,sigma_ext_forcing)
    if update_mode in [:roa,:ga]
        error("Not yet implemented for random asyncrhonous updating")
    end
    while i < nrand_init
        rand_init_net!(net;update_mode=update_mode,init_h=:rand,init_constraints...)
        p_att = Vector{Tuple{Vararg{Tuple{Vararg{Int8,net_size}},net_hist_size}}}()
        init_probe = tuple(map(x->tuple(x...),(net.sim_buffer[:,ti] for ti in 1:net_hist_size))...)
        push!(p_att,init_probe)
        add_attractor = false
        already_found = false
        attractor_found = false
        tmax = 0
        for t in 0:max_steps-1
            hX[:] = evolve_net!(net,1,keep=false,update_mode=update_mode,sim=sim_att_buffer,t_init=t,forcing_rhythms=forcing_rhythms)
            ntp1 = tuple(map(x->tuple(x...),(hX[:,ti] for ti in 1:net_hist_size))...)
            tmax = t
            for nti in collect(size(p_att,1):-1:1)
                if ntp1[end] == p_att[nti][end]
                    #println("Likely attractor found:$i at t=$t")
                    # Check if it was previously found
                    for (pi,att) in enumerate(att_list)
                        for p_att_i in p_att
                            if p_att_i in att
                                already_found = true
                                if !(init_probe in att)
                                    basin_sizes[pi]+=1
                                    push!(transients[pi],nti-1)
                                end
                                break
                            end
                        end
                        if already_found
                            break
                        end
                    end
                    if !already_found
                        deleteat!(p_att,1:nti-1)
                        push!(att_list,p_att)
                        push!(basin_sizes,init_probe in p_att?0:1)
                        push!(transients,[nti-1])
                        println("Attractor ",size(att_list,1)," has period ",size(p_att,1))
                        for att_i in p_att
                            println(join(att_i[end],","))
                        end
                    end
                    attractor_found = true
                    break
                end
            end
            if attractor_found
                #println(i," : ",nti-1)
                break
            end
            push!(p_att,ntp1)
        end
        i += 1
        if !(already_found||attractor_found)
            println("No attractor found After $tmax steps")
            print("probe:")
            println(join(init_probe[end],","))
            #  println("No attractors after reaching maxtimum time steps")
            unsuccess += 1
            # println("No attractor found for initial condition $i:")
            # for probe_i in probe
            #     conf = reverse(digits(probe_i,tag_base,net_size))
            #     println(join(conf,","))
            # end
        end
    end
    if length(att_list)!=0
        println("Estimated basin sizes:")
        for (pi,att) in enumerate(att_list)
            println("Attractor $pi:",format(100.0*basin_sizes[pi]/sum(basin_sizes),precision=2,zeropadding=true),"%")
        end
    end
    println(unsuccess," unsuccessful tries after reaching maximmum time steps set")
    return att_list, transients
end
