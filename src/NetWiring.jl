using LightGraphs
using IterTools

function load_net2(dir::AbstractString)
    nodeFfiles = readdir(dir)
    filter!(x->ismatch(r"\w+.txt$",x),nodeFfiles)
    node_ids_byname = Dict{AbstractString,Integer}()
    node_ids_bynum = Dict{Integer,AbstractString}()
    node_tables = Dict{AbstractString,Array{AbstractString,1}}()

    for name in	nodeFfiles
        f = open(joinpath(dir,name))
        try
            node_id = split(name,['.','_'])[1:2]
            node_tables[node_id[2]] = readlines(f)
            node_number = parse(Int,node_id[1])
            node_ids_byname[node_id[2]] = node_number
            node_ids_bynum[node_number] = node_id[2]
        finally
            close(f)
        end
    end

    netstate = Vector{BNode}()
    for (node_num,node_name) in sort(collect(node_ids_bynum))
        push!(netstate,BNode(node_name))
    end
    net = Net(netstate,DiGraph(length(netstate)))
    for sigma in net.state
        test_pattern = match(r"\d+", strip(node_tables[sigma.id][1]))
        sigma.def_state = test_pattern != nothing ? parse(Int,test_pattern.match) : error("Invalid value: First line of node file $(sigma.id) must contain a numeric value to set default state")
        # set_node!(sigma)
        test_pattern = match(r"\d+", strip(node_tables[sigma.id][2]))
        max_states = test_pattern != nothing ? parse(Int,test_pattern.match) : error("Invalid value: Second line of node file $(sigma.id) must contain the maximum number of possible states")
        sigma.state_range = collect(range(0, max_states)) # All possible states is given by a sequence of consecutive integer numbers starting from 0
        push!(sigma.thresholds,(max_states-1)/max_states*sigma.state_range[2:end]...)
        reg_ids = map(x->split(x,['_','.']),matchall(r"(\d+_)?[[:alnum:]]+(\.\d+)?",strip(node_tables[sigma.id][3])))
        for reg_id in reg_ids
            if length(reg_id) == 1
                reg_id_name = reg_id[1]
            elseif length(reg_id) == 2
                reg_id_name = reg_id[2]
                push!(sigma.tau_in,0)
            elseif length(reg_id) == 3
                reg_id_name = reg_id[2]
                push!(sigma.tau_in,parse(Int,reg_id[3]))
            else
                error("Format error for list of regulators for node $(sigma.id)")
            end
            index_reg = findfirst(x -> x.id==reg_id_name,net.state)
            index_reg!=0 || error("Unregistered node in list of regulators: There is no file for node $reg_id_name")
            push!(sigma.sigma_in, index_reg)
            add_edge!(net.graph,findfirst(x -> x.id==reg_id_name,net.state),findfirst(x -> x.id==sigma.id,net.state))
        end
        for reg_phrase in node_tables[sigma.id][4:end]
            r = matchall(r"\d+",reg_phrase)
            if(!isempty(r))
                sigma.F[tuple([parse(Int,ri) for ri in r[1:(end-1)]]...)] = parse(Int,r[end])
            end
        end
    end
    for n in vertices(net.graph)
        push!(net.state[n].sigma_out,out_neighbors(net.graph,n)...)
    end
    mem_size = maximum(vcat([sigma.tau_in for sigma in net.state]...)) + 1
    for (nsigma,sigma) in enumerate(net.state)
        # sigma.state = fill(sigma.def_state,maximum([sigmak.tau_in[findfirst(sigmak.sigma_in,nsigma)] for sigmak in net.state[sigma.sigma_out]])+1)
        sigma.state = fill(sigma.def_state,mem_size)
        sigma.sigma_in_state = [net.state[sigmak].def_state for sigmak in sigma.sigma_in]
        sigma.sigma_in_state_mu = Vector{Float16}(length(sigma.sigma_in))
        # set_node!(sigma)
    end
    return net
end

function fill_tables!(netstate::Vector{BNode},p_scheme::Symbol,p::AbstractFloat)
    if p_scheme==:semifix
        for sigma in netstate
            itp = product((netstate[sigmak].state_range for sigmak in sigma.sigma_in)...)
            table_size = length(itp)
            n_ones = round(Int,p*table_size)
            table_output = trues(n_ones)
            push!(table_output,bitrand(table_size-n_ones)...)
            shuffle!(table_output)
            sigma.F = Dict(reg_phrase => table_output[i] for (i,reg_phrase) in enumerate(itp))
        end
    elseif p_scheme==:prob
        for sigma in netstate
            itp = product((netstate[sigmak].state_range for sigmak in sigma.sigma_in)...)
            table_output = rand(length(itp)) .< p
            sigma.F = Dict(reg_phrase => table_output[i] for (i,reg_phrase) in enumerate(itp))
        end
    elseif p_scheme==:fix
        for sigma in netstate
            itp = product((netstate[sigmak].state_range for sigmak in sigma.sigma_in)...)
            table_size = length(itp)
            table_output = falses(table_size)
            n_ones = round(Int,p*table_size)
            table_output[1:n_ones] = true
            shuffle!(table_output)
            sigma.F = Dict(reg_phrase => table_output[i] for (i,reg_phrase) in enumerate(itp))
        end
    else
        error("Invalid scheme name: Choose from :fix,:prob,:semifix")
    end
    return nothing
end

function fill_taus!(net_graph::DiGraph,netstate::Vector{BNode},tau_scheme::Symbol;tau_params...)
    tau = 0
    taus = [0]
    found = false
    if tau_scheme==:fix
        for (id,value) in tau_params
            if id==:tau
                tau = value
                found = true
                break
            end
        end
        found || error("Expected tau parameter")
        for (nsigma,sigma) in enumerate(netstate)
            push!(sigma.tau_in,fill(tau,length(sigma.sigma_in))...)
            push!(sigma.sigma_in_state,zeros(Int8,length(sigma.sigma_in))...)
            sigma.sigma_in_state_mu = Vector{Float16}(length(sigma.sigma_in))
            # set_node!(sigma)
        end
    elseif tau_scheme==:random
        for (id,value) in tau_params
            if id==:taus
                taus = value
                found = true
                break
            end
        end
        found || error("Expected taus parameter")
        for (nsigma,sigma) in enumerate(netstate)
            push!(sigma.tau_in,rand(taus,length(sigma.sigma_in))...)
            push!(sigma.sigma_in_state,zeros(Int8,length(sigma.sigma_in))...)
            sigma.sigma_in_state_mu = Vector{Float16}(length(sigma.sigma_in))
            # set_node!(sigma)
        end
    else
        error("Invalid scheme name: Choose from :fix,:random")
    end
    mem_size = maximum(vcat([sigma.tau_in for sigma in netstate]...)) + 1
    for sigma in netstate
        sigma.state = zeros(Int8,mem_size)
    end
    return nothing
end

function generate_RBN(N::Integer,K::Integer;p_scheme::Symbol=:fix,tau_scheme::Symbol=:fix,params...)
    p = 0.5
    tau = 0
    taus = [0]
    for (id,value) in params
        if id==:p
            p = value
        elseif id==:tau
            tau = value
        elseif id==:taus
            taus = value
        end
    end
    net_graph = random_regular_digraph(N, K; dir=:in, seed=-1)
    netstate = Vector{BNode}(N)
    findfirst([:fix,:prob,:semifix],p_scheme)!=0 || error("Invalid scheme name: Choose from :fix,:prob,:semifix")
    for n in vertices(net_graph)
        netstate[n] = BNode(join(["n",string(n)]))
    end
    println("Succesfully added ",length(netstate)," nodes")
    for n in vertices(net_graph)
        isempty(in_neighbors(net_graph,n))||push!(netstate[n].sigma_in, in_neighbors(net_graph,n)...)
        isempty(out_neighbors(net_graph,n))||push!(netstate[n].sigma_out, out_neighbors(net_graph,n)...)
    end
    fill_tables!(netstate,p_scheme,p)
    fill_taus!(net_graph,netstate,tau_scheme,tau=tau,taus=taus)
    return Net(netstate,net_graph)
end

function knock_out!(netstate::Vector{BNode},nodeKO_id::AbstractString)
    index_KO = findfirst(x->x.id==nodeKO_id,netstate)
    if index_KO!=0
        nodeKO = netstate[index_KO]
        init_node!(nodeKO,nodeKO.def_state)
        nodeKO.thresholds = [nodeKO.thresholds[findfirst(sort([nodeKO.def_state,nodeKO.thresholds...]),nodeKO.def_state)]]
        nodeKO.state_range = [nodeKO.def_state]
        for k in keys(nodeKO.F)
            nodeKO.F[k] = nodeKO.def_state
        end
    else
        error("Node id invalid: $nodeKO_id does not exist in the network")
    end
end

function knock_down!(netstate::Vector{BNode},nodeKD_id::AbstractString)
    index_KD = findfirst(x->x.id==nodeKD_id,netstate)
    if index_KD!=0
        nodeKD = netstate[index_KD]
        init_node!(nodeKD,nodeKO.def_state)
        for k in keys(nodeKD.F)
            nodeKD.F[k] = nodeKD.def_state
        end
    else
        error("Node id invalid: $nodeKD_id does not exist in the network")
    end
end

function constitutive_state!(netstate::Vector{BNode},const_node_id::AbstractString,const_state::Integer)
    const_node = filter(x->x.id==const_node_id,netstate)
    if(!isempty(const_node))
        const_node = const_node[1]
        [const_node.F[k] = const_state for k in keys(const_node.F)]
        init_node!(const_node,const_state)
        #const_node.state_range = [const_state]
    else
        error("Node id invalid: $nodeCA_id does not exist in the network")
    end
end
