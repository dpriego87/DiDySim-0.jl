using LightGraphs
using IterTools

function load_net(dir::AbstractString)
    nodeFfiles = readdir(dir)
    filter!(x->ismatch(r"\w+.txt$",x),nodeFfiles)
    # node_ids_byname = Dict{AbstractString,Integer}()
    node_ids_bynum = Dict{Integer,AbstractString}()
    node_tables = Dict{AbstractString,Array{AbstractString,1}}()

    for name in	nodeFfiles
        f = open(joinpath(dir,name))
        try
            node_id = split(name,['.','_'])[1:2]
            node_tables[node_id[2]] = readlines(f)
            node_number = parse(Int,node_id[1])
            # node_ids_byname[node_id[2]] = node_number
            node_ids_bynum[node_number] = String(node_id[2])
        finally
            close(f)
        end
    end
    allunique(values(node_ids_bynum))||error("Node identifiers must be unique")
    ids = [node_id for (node_num,node_id) in sort(collect(node_ids_bynum))]
    nodes = Vector{BNode2{Int8}}(length(node_ids_bynum))
    graph = DiGraph(length(node_ids_bynum))
    tau_matrix::SparseMatrixCSC{Int,Int} = spzeros(nv(graph),nv(graph))
    for (nsigma,id) in enumerate(ids)
        test_pattern = matchall(r"\d+", strip(node_tables[id][1]))
        !isempty(test_pattern)||error("Invalid value: First line of node file $id must contain a numeric value to set default state (optional extra field, separated by space, order number for semi-synchronous updating)")
        if length(test_pattern)==1
            def_state,sync_order = parse(Int,test_pattern[1]),1
        elseif length(test_pattern)==2
            def_state,sync_order = map(x->parse(Int,x),test_pattern)
        else
            error("Incorrect number of arguments: Less or equal than 2 arguments expected")
        end
        test_pattern = match(r"\d+", strip(node_tables[id][2]))
        max_states = test_pattern != nothing ? parse(Int,test_pattern.match) : error("Invalid value: Second line of node file $id must contain the maximum number of possible states")
        state_range = collect(range(0, max_states)) # All possible states is given by a sequence of consecutive integer numbers starting from 0
        sigma = nodes[nsigma] = BNode2{Int8}(id,def_state,state_range,sync_order)
        push!(sigma.thresholds,(max_states-1)/max_states*sigma.state_range[2:end]...)
        reg_ids = map(x->split(x,['_','.']),matchall(r"(\d+_)?[[:alnum:]]+(\.\d+)?",strip(node_tables[id][3])))
        reg_indices = Int[]
        for reg_id in reg_ids
            if length(reg_id) == 1 # If only the id is typed,
                reg_id_name = reg_id[1]
            elseif length(reg_id) == 2
                if ismatch(r"^\d+$",reg_id[1])
                    reg_id_name = reg_id[2]
                    tau = 0
                else
                    reg_id_name = reg_id[1]
                    tau = parse(reg_id[2])
                end
            elseif length(reg_id) == 3
                reg_id_name = reg_id[2]
                tau = parse(Int,reg_id[3])
            else
                error("Format error for list of regulators for node ",sigma.id)
            end
            index_reg = findfirst(ids,String(reg_id_name))
            index_reg!=0 || error("Unregistered node in list of regulators: There is no file for node $reg_id_name")
            push!(reg_indices,index_reg)
            tau_matrix[index_reg,nsigma] = tau
            add_edge!(graph,index_reg,nsigma)
        end
        sorted_reg_indices = sortperm(reg_indices)
        for reg_phrase in node_tables[sigma.id][4:end]
            r = matchall(r"\d+",reg_phrase)
            if !isempty(r)
                sigma.rule[[parse(Int,ri) for ri in r[1:(end-1)]][sorted_reg_indices]] = parse(Int,r[end])
            end
        end
    end
    if isempty(nonzeros(tau_matrix))
        return Net2{Int8}(nodes,graph)
    else
        return Net2{Int8}(nodes,graph,tau_matrix)
    end
end

function to_boolean!(net::Net2{Int8})
    netstateB = Vector{BNode2{Bool}}()
    ids = [sigma.id for sigma in net.nodes]
    #Add auxiliary nodes
    for (ni,sigma) in enumerate(net.nodes)
        max_state = length(sigma.state_range)
        if max_state>2
            def_state_index = findfirst(sort(sigma.state_range),sigma.def_state)
            max_digits = length(digits(max_state,2)) #transform max_state to base two to determine maximum nodes needed to reencode multistate nodes
            for s in max_digits:-1:1
                push!(netstateB,BNode2{Bool}(join([sigma.id,s-1,s]),digits(def_state_index-1,2,max_digits)[s],[false,true],sigma.sync_order))
                push!(sigma.thresholds,0.5) #Watch out, this might result in a big difference with respect to the original multistate node
            end
        else
            push!(netstate,BNode2{Bool}(sigma.id,sigma.def_state==1,[false,true],sigma.sync_order))
        end
    end
    graphB = DiGraph(length(netstateB))
    #Rewire network
    idsB = [sigma.id for sigma in netstateB]
    indexB = 0
    for (ni,sigma) in enumerate(net.nodes)
        max_state = length(sigma.state_range)
        max_digits = length(digits(max_state,2))
        in_neighbors(net.graph,ni)

        for reg_id in reg_ids
            if length(reg_id) == 1 # If only the id is typed,
                reg_id_name = reg_id[1]
            elseif length(reg_id) == 2
                if ismatch(r"^\d+$",reg_id[1])
                    reg_id_name = reg_id[2]
                    tau = 0
                else
                    reg_id_name = reg_id[1]
                    tau = parse(reg_id[2])
                end
            elseif length(reg_id) == 3
                reg_id_name = reg_id[2]
                tau = parse(Int,reg_id[3])
            else
                error("Format error for list of regulators for node ",sigma.id)
            end
            index_reg = findfirst(ids,reg_id_name)
            index_reg!=0 || error("Unregistered node in list of regulators: There is no file for node $reg_id_name")
            push!(reg_indices,index_reg)
            tau_matrix[index_reg,nsigma] = tau
            add_edge!(graph,index_reg,nsigma)
        end
        sorted_reg_indices = sortperm(reg_indices)
        for reg_phrase in node_tables[sigma.id][4:end]
            r = matchall(r"\d+",reg_phrase)
            if !isempty(r)
                sigma.rule[[parse(Int,ri) for ri in r[1:(end-1)]][sorted_reg_indices]] = parse(Int,r[end])
            end
        end
        indexB = indexB + max_digits
    end

    netB = Net2{Bool}(netstate,graph)

    return netB
end


function generate_RBN2(N::Integer,K::Integer;p_scheme::Symbol=:fix,tau_scheme::Symbol=:fix,params...)
    findfirst([:fix,:prob,:semifix],p_scheme)!=0 || error("Invalid scheme name: Choose from :fix,:prob,:semifix")
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
    net_nodes = Vector{BNode2{Bool}}(N)
    for n in vertices(net_graph)
        net_nodes[n] = BNode2{Bool}(join(["n",string(n)]),rand([false,true]),[false,true])
    end
    println("Succesfully created ",length(net_nodes)," nodes")
    fill_tables!(net_nodes,net_graph,p_scheme,p)
    tau_matrix = fill_taus(net_graph,tau_scheme,tau=tau,taus=taus)
    return Net2{Bool}(net_nodes,net_graph,tau_matrix)
end

function fill_tables!(net_nodes::Vector{BNode2{Bool}},net_graph::DiGraph,p_scheme::Symbol,p::AbstractFloat)
    A = adjacency_matrix(net_graph)
    rows = rowvals(A)
    if p_scheme==:semifix
        for (ns,sigma) in enumerate(netnet_nodes)
            itp = product((net_nodes[sigmak].state_range for sigmak in rows[nzrange(A,ns)])...)
            table_size = length(itp)
            n_ones = round(Int,p*table_size)
            table_output = trues(n_ones)
            push!(table_output,bitrand(table_size-n_ones)...)
            shuffle!(table_output)
            sigma.rule = Dict(BitArray(collect(reg_phrase)) => table_output[i] for (i,reg_phrase) in enumerate(itp))::reg_function
        end
    elseif p_scheme==:prob
        for (ns,sigma) in enumerate(net_nodes)
            itp = product((net_nodes[sigmak].state_range for sigmak in rows[nzrange(A,ns)])...)
            table_output = rand(length(itp)) .< p
            sigma.rule = Dict(BitArray(collect(reg_phrase)) => table_output[i] for (i,reg_phrase) in enumerate(itp))::reg_function
        end
    elseif p_scheme==:fix
        for (ns,sigma) in enumerate(net_nodes)
            itp = product((net_nodes[sigmak].state_range for sigmak in rows[nzrange(A,ns)])...)
            table_size = length(itp)
            table_output = falses(table_size)
            n_ones = round(Int,p*table_size)
            table_output[1:n_ones] = true
            shuffle!(table_output)
            sigma.rule = Dict(BitArray(collect(reg_phrase)) => table_output[i] for (i,reg_phrase) in enumerate(itp))::reg_function
        end
    else
        error("Invalid scheme name: Choose from :fix,:prob,:semifix")
    end
    return nothing
end

function fill_taus(net_graph::DiGraph,tau_scheme::Symbol;tau_params...)
    tau = 0
    taus = [0]
    found = false
    tau_matrix::SparseMatrixCSC{Int,Int} = spzeros(nv(net_graph),nv(net_graph))
    if tau_scheme==:fix
        for (id,value) in tau_params
            if id==:tau
                tau = value
                found = true
                break
            end
        end
        found || error("Expected tau parameter")
        if tau != 0
            tau_matrix = adjacency_matrix(net_graph)
            tau_vals = nonzeros(tau_matrix)
            tau_vals[:] = tau
        end
        # for (nsigma,sigma) in enumerate(netstate)
        #     push!(sigma.tau_in,fill(tau,length(sigma.sigma_in))...)
        #     push!(sigma.sigma_in_state,zeros(Int8,length(sigma.sigma_in))...)
        #     sigma.sigma_in_state_mu = Vector{Float16}(length(sigma.sigma_in))
        # end
    elseif tau_scheme==:random
        for (id,value) in tau_params
            if id==:taus
                taus = value
                found = true
                break
            end
        end
        found || error("Expected taus parameter")
        tau_matrix = adjacency_matrix(net_graph)
        tau_vals = nonzeros(tau_matrix)
        tau_vals[:] = rand(taus,length(tau_vals))
        # for (nsigma,sigma) in enumerate(netstate)
        #     push!(sigma.tau_in,rand(taus,length(sigma.sigma_in))...)
        #     push!(sigma.sigma_in_state,zeros(Int8,length(sigma.sigma_in))...)
        #     sigma.sigma_in_state_mu = Vector{Float16}(length(sigma.sigma_in))
        # end
    else
        error("Invalid scheme name: Choose from :fix,:random")
    end
    # mem_size = maximum(vcat([sigma.tau_in for sigma in netstate]...)) + 1
    # for sigma in netstate
    #     sigma.state = zeros(Int8,mem_size)
    # end
    return tau_matrix
end

function knock_out!(net::Net2,nodeKO_id::Union{AbstractString,Int})
    index_KO = typeof(nodeKO_id)<:Int ? nodeKO_id : findfirst(sigma->sigma.id==nodeKO_id,net.nodes)
	0 < index_KO <= length(net.nodes) || error("Unexisting node: $nodeKO_id")
    nodeKO = net.nodes[index_KO]
    net.sim_buffer[index_KO,:] = nodeKO.def_state
    # nodeKO.thresholds = [nodeKO.thresholds[findfirst(sort([nodeKO.def_state,nodeKO.thresholds...]),nodeKO.def_state)]]
    nodeKO.state_range = [nodeKO.def_state]
    for k in keys(nodeKO.rule)
        nodeKO.rule[k] = nodeKO.def_state
    end
end

function knock_down!(net::Net2,nodeKD_id::Union{AbstractString,Int})
    index_KD = typeof(nodeKD_id)<:Int ? nodeKD_id : findfirst(sigma->sigma.id==nodeKD_id,net.nodes)
	0 < index_KD <= length(net.nodes) || error("Unexisting node: $nodeKD_id")
    nodeKD = net.nodes[index_KD]
    net.sim_buffer[index_KD,:] = nodeKD.def_state
    for k in keys(nodeKD.rule)
        nodeKD.rule[k] = nodeKD.def_state
    end
end

function fix_state!(net::Net2,const_node_id::Union{AbstractString,Int},const_state::Integer)
    index_const = typeof(const_node_id)<:Int ? const_node_id : findfirst(sigma->sigma.id==const_node_id,net.nodes)
	0 < index_const <= length(net.nodes) || error("Unexisting node: $const_node_id")
    node_const = net.nodes[index_const]
    const_state in node_const.state_range || error("Trying to set node $(node_const.id) to a forbidden state")
    net.sim_buffer[index_const,:] = const_state
    for k in keys(node_const.rule)
        node_const.rule[k] = const_state
    end
end
