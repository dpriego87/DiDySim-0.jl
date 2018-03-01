using StatsBase
using GZip

macro update_schemes()
	:(Dict(:sy=>:sy,:fd=>fixed_delay!,:md=>smoothed_fixed_delay!,:st=>strob!,:mst=>smoothed_strob!,:ga=>:ga,:roa=>:roa,:ssy=>semi_sync!))
end
function init_net!(net::Net2,node_indices::Vector{Int}=collect(vertices(net.graph));update_mode::Symbol=:sy,t0_mode::Symbol=:custom,init_h::Symbol=:last,init_array::Vector{Int}=Int[],init_constraints...)
	allunique(node_indices) || error("Duplicate node indices")
	intersect(vertices(net.graph),node_indices)==node_indices || error("Invalid node indices")
	umodes = @update_schemes
	haskey(umodes,update_mode) || error("Unexisting update scheme. Available schemes are :",keys(umodes))
	nodes_to_constraint = Int[]
	max_buffer_tlength = size(net.sim_buffer)[2]
	if !isempty(init_constraints)
		ids = [string(id) for (id,state) in init_constraints]
		allunique(ids) || error("Duplicate constraint specification")
		intersect(ids,(sigma.id for sigma in net.nodes))==ids || error("Trying to set unexisting nodes:",ids)
		all(i->i<:Union{Int,Bool,Vector{Int},Vector{Bool}},(typeof(state) for (id,state) in init_constraints)) || error("Wrong format to specify state")
		push!(nodes_to_constraint,(findfirst([sigma.id for sigma in net.nodes],id) for id in ids)...)
		all(map((x,y)->typeof(x)<:Vector?intersect(x,y)==x : intersect([x],y)==[x],(state for (id,state) in init_constraints),(sigma.state_range for sigma in net.nodes[nodes_to_constraint]))) || error("Trying to set node to a forbidden state")
		node_indices = setdiff(node_indices,nodes_to_constraint)
		states = Vector{Vector{Int}}(length(init_constraints))
		for (ni,nsigma) in enumerate(nodes_to_constraint)
			state = init_constraints[ni][2]
			states[ni] = isempty(state) ? [net.nodes[nsigma].def_state] : (typeof(state)<:Vector ? state : [state])
		end
	end
	if max_buffer_tlength > 1
		if init_h==:default
			for ni in vcat(node_indices,nodes_to_constraint)
				net.sim_buffer[ni,1:end-1] = net.nodes[ni].def_state
			end
		elseif init_h==:rand
			for ni in node_indices
				net.sim_buffer[ni,1:end-1] = rand(net.nodes[ni].state_range,max_buffer_tlength-1)
			end
		elseif init_h!=:last
			error("Invalid initialization flag. Available options are\n:default for default, :rand for random")
		end
	end
	if t0_mode==:default
		net.sim_buffer[node_indices,end] = [sigma.def_state for sigma in net.nodes[node_indices]]
	elseif t0_mode==:custom
		if !isempty(init_array)
			length(init_array)==length(node_indices)||error("Initial configuration vector must have the same length of nodes vector")
			for (n,sigma) in enumerate(net.nodes[node_indices])
				if init_array[n] in sigma.state_range
					net.sim_buffer[node_indices[n],end] = init_array[n]
				else
					error("Trying to set node ",sigma.id," to a forbidden state, value ",init_array[n]," is out of possible states range ",sigma.state_range)
				end
			end
		end
	else
		error("Invalid initialization flag. Available options are\n:default for default, :rand for random :custom start from last available configuration")
	end
	for (ni,nsigma) in enumerate(nodes_to_constraint)
		state = states[ni]
		if length(state)<=max_buffer_tlength
			net.sim_buffer[ni,end-length(state)+1:end] = collect(state)
		else
			net.sim_buffer[ni,:] = collect(state[end-max_buffer_tlength+1:end])
		end
	end

	in_state_vals = nonzeros(net.in_states)
	# mu_in_state_vals = nonzeros(net.mu_in_states)
	nsigma_ins = rowvals(net.in_states)
	# tau_vals = nonzeros(net.in_taus)
	# if update_mode == :sy
	# 	for i in vertices(net.graph)
	# 		# j_in = nzrange(net.in_states,i)
	# 		# mu_in_state_vals[j_in] = in_state_vals[j_in] = net.sim_buffer[nsigma_ins[j_in],end]
	# 	end
	if update_mode == :fd || update_mode == :st
		for i in vertices(net.graph)
			for j in nzrange(net.in_states,i)
				in_state_vals[j] = net.sim_buffer[nsigma_ins[j],end-net.in_taus[nsigma_ins[j],i]]
			end
		end
	elseif update_mode == :md || update_mode == :mst
		mu_in_state_vals = nonzeros(net.mu_in_states)
		for i in vertices(net.graph)
			for j in nzrange(net.in_states,i)
				nsigma = nsigma_ins[j]
				tau = net.in_taus[nsigma,i]
				sigma_in = net.nodes[nsigma]
				if tau==0
					mu_in_state_vals[j] = in_state_vals[j] = net.sim_buffer[nsigma,end]
				else
					mu_in_state_vals[j] = mean_sigmak = mean(net.sim_buffer[nsigma,end-tau:end])
					in_state_vals[j] = sigma_in.state_range[findlast(mean_sigmak .>= sigma_in.thresholds)+1]
				end
			end
		end
	end
	return net.sim_buffer[:,end]
end

function rand_init_net!(net::Net2,node_indices::Vector{Int}=collect(vertices(net.graph));update_mode::Symbol=:sy,init_h::Symbol=:rand,init_constraints...)
	allunique(node_indices) || error("Duplicate node indices")
	intersect(vertices(net.graph),node_indices)==node_indices || error("Invalid node indices")
	umodes = @update_schemes
	haskey(umodes,update_mode) || error("Unexisting update scheme. Available schemes are :",keys(umodes))
	nodes_to_constraint = Int[]
	max_buffer_tlength = size(net.sim_buffer)[2]
	if !isempty(init_constraints)
		ids = [string(id) for (id,state) in init_constraints]
		allunique(ids) || error("Duplicate constraint specification")
		intersect(ids,(sigma.id for sigma in net.nodes))==ids || error("Trying to set unexisting node $(setdiff(ids,(sigma.id for sigma in net.nodes)))")
		all(i->i<:Union{Int,Bool,Vector{Int},Vector{Bool}},(typeof(state) for (id,state) in init_constraints)) || error("Wrong format to specify state")
		push!(nodes_to_constraint,(findfirst([sigma.id for sigma in net.nodes],id) for id in ids)...)
		all(map((x,y)->typeof(x)<:Vector?intersect(x,y)==x : intersect([x],y)==[x],(state for (id,state) in init_constraints),(sigma.state_range for sigma in net.nodes[nodes_to_constraint]))) || error("Trying to set node to a forbidden state")
		node_indices = setdiff(node_indices,nodes_to_constraint)
		states = Vector{Vector{Int}}(length(init_constraints))
		for (ni,nsigma) in enumerate(nodes_to_constraint)
			state = init_constraints[ni][2]
			states[ni] = isempty(state) ? [net.nodes[nsigma].def_state] : unique(state)
			net.sim_buffer[nsigma,end] = length(state)>1 ? rand(state) : state[1]
		end
	end
	net.sim_buffer[node_indices,end] = [rand(sigma.state_range) for sigma in net.nodes[node_indices]]
	if max_buffer_tlength > 1
		if init_h==:default
			for ni in vcat(node_indices,nodes_to_constraint)
				net.sim_buffer[ni,1:end-1] = net.nodes[ni].def_state
			end
		elseif init_h==:rand
			for ni in node_indices
				net.sim_buffer[ni,1:end-1] = rand(net.nodes[ni].state_range,max_buffer_tlength-1)
			end
			for (n,ni) in enumerate(nodes_to_constraint)
				net.sim_buffer[ni,1:end-1] = rand(states[n],max_buffer_tlength-1)
			end
		elseif init_h!=:last
			error("Invalid initialization flag. Available options are\n:default for default, :rand for random")
		end
	end
	in_state_vals = nonzeros(net.in_states)
	nsigma_ins = rowvals(net.in_states)
	# tau_vals = nonzeros(net.in_taus)
	if update_mode == :fd || update_mode == :st
		for i in vertices(net.graph)
			for j in nzrange(net.in_states,i)
				in_state_vals[j] = net.sim_buffer[nsigma_ins[j],end-net.in_taus[nsigma_ins[j],i]]
			end
		end
	elseif update_mode == :md || update_mode == :mst
		mu_in_state_vals = nonzeros(net.mu_in_states)
		for i in vertices(net.graph)
			for j in nzrange(net.in_states,i)
				nsigma = nsigma_ins[j]
				tau = net.in_taus[nsigma,i]
				sigma_in = net.nodes[nsigma]
				if tau==0
					mu_in_state_vals[j] = in_state_vals[j] = net.sim_buffer[nsigma,end]
				else
					mu_in_state_vals[j] = mean_sigmak = mean(net.sim_buffer[nsigma,end-tau:end])
					in_state_vals[j] = sigma_in.state_range[findlast(mean_sigmak .>= sigma_in.thresholds)+1]
				end
			end
		end
	end
	# return net.sim_buffer[:,end]
	return nothing
end

macro updating_scheme(x)
	update_ex0 = quote
	end
	if eval(x) == :sy
		update_ex1 = quote
			setdiff(sigma_indices,forced_indices)
		end
		update_ex2 = quote
			sim[i,t+1] = net.nodes[i].rule[sim[nsigma_ins[nzrange(net.in_states,i)],t]]
		end
	elseif eval(x) == :ga
		update_ex1 = quote
			rand(sigma_indices)
		end
		update_ex2 = quote
			unchanged = setdiff(sigma_indices,forced_indices)
			if !(i in forced_indices)
				sim[i,t+1] = net.nodes[i].rule[sim[nsigma_ins[nzrange(net.in_states,i)],t]]
				unchanged = setdiff(unchanged,i)
			end
			# @show unchanged
			sim[unchanged,t+1] = sim[unchanged,t]
		end
	elseif eval(x) ==:roa
		update_ex0 = quote
			sim[:,t+1] = sim[:,t]
		end
		update_ex1 = quote
			shuffle(sigma_indices)
		end
		update_ex2 = quote
			if !(i in forced_indices)
				sim[i,t+1] = net.nodes[i].rule[sim[nsigma_ins[nzrange(net.in_states,i)],t+1]]
			end
		end
	end
	return esc(quote
		nsigma_ins = rowvals(net.in_states)
		sigma_noise_indices = [ni for (ni,eta) in sigma_noises]
		sigma_indices = 1:length(net.nodes)
		for t in i_t0:i_t0+steps-1
			$update_ex0
			forced_indices = Int[]
			# Updating forced nodes
			for (nsigma,rhythm_values) in nsigma_ext_forcing
				rhythm,forced_values = rhythm_values
				index_rhythm = mod(t_init+t-i_t0,length(rhythm))+1
				if rhythm[index_rhythm]
					#println("At t=",t_init+t+1-i_t0,",trying to fix node index ",nsigma," on rhythm index ",index_rhythm," to forced_value ",forced_values[1])
					sim[nsigma,t+1] = forced_values[1]
					push!(forced_indices,nsigma)
					nsigma_ext_forcing[nsigma][2] = circshift(forced_values,-1)
				end
			end
			# Apply the original rule except for the forced nodes
			for i in $update_ex1
				$update_ex2
				# if i in forced_indices
				# 	println("At t=",t_init+t+1-i_t0,",node index ",i," changed from ",sim[i,t]," to ",sim[i,t+1])
				# end
				# Add noise where specified, expect for those forced nodes
				if i in sigma_noise_indices
					eta = sigma_noises[findfirst(sigma_noise_indices,i)][2]
					true_value = sim[i,t+1]
					sim[i,t+1] = rand()<=eta ? rand(filter(x -> x!=true_value,net.nodes[i].state_range)) : true_value
				end
				# sim[i,t+1] = net.nodes[i].rule[sim[nsigma_ins[nzrange(net.in_states,i)],t]]
			end
		end
	end)
end

function evolve_net!{T<:state_type}(net::Net2{T},steps::Int;keep::Bool=true,update_mode::Symbol=:sy,sim::Array{T,2}=Array{T}(nv(net.graph),steps + 1 + (update_mode == :sy ? 0 : maximum(net.in_taus))),t_init::Int=0,noise_vector::Vector{Pair{Symbol,Float64}}=Pair{Symbol,Float64}[],forcing_rhythms::Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}=Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}())
	umodes = @update_schemes
	haskey(umodes,update_mode) || error("Unexisting update scheme. Available schemes are :",keys(umodes))
	steps>=0 || error("Steps must be a positive integer")
	hist_size = 1 + (update_mode == :sy ? 0 : maximum(net.in_taus))
	# Initialize simulation table (from t_init - tau0 to t_init)
	size(sim) == (nv(net.graph),steps + hist_size) || error("Array for storing results has incorrect dimensions ",size(sim), ". Expected array size is ",(nv(net.graph),steps + hist_size))
	if update_mode == :sy
		sim[:,1] = net.sim_buffer[:,end]
	else
		sim[:,1:hist_size] = net.sim_buffer[:]
	end
	sigma_noises = collect(findfirst(x->x.id==string(sigma),net.nodes)=>eta for (sigma,eta) in noise_vector)
	nsigma_ext_forcing = Dict(findfirst(x->x.id==string(sigma_symbol),net.nodes)=>collect(rhythm_states) for (sigma_symbol,rhythm_states) in forcing_rhythms)
	i_t0 = hist_size
	if update_mode == :sy
		@updating_scheme :sy
	elseif update_mode == :ga
		for (nsigma,rhythm_states) in nsigma_ext_forcing
			nsigma_ext_forcing[nsigma][1] = repeat(nsigma_ext_forcing[nsigma][1],inner=[length(net.nodes)])
			nsigma_ext_forcing[nsigma][2] = repeat(nsigma_ext_forcing[nsigma][2],inner=[length(net.nodes)])
		end
		@show nsigma_ext_forcing
		@updating_scheme :ga
	elseif update_mode == :roa
		@updating_scheme :roa
	else
		umodes[update_mode](net,sim,steps,hist_size,sigma_noises,sigma_ext_forcing,t_init)
	end
	buffer_size = size(net.sim_buffer)[2]
	if buffer_size == hist_size || steps > buffer_size
		net.sim_buffer[:] = sim[:,end-buffer_size+1:end]
	elseif steps <= hist_size
		net.sim_buffer[:,1:end-steps-hist_size] = net.sim_buffer[:,1+steps+hist_size:end]
		net.sim_buffer[:,end-steps-hist_size+1:end] = sim[:]
	else
		@show steps, hist_size,buffer_size
		error("Wrong dimensions")
	end
	if keep
		return sim
	else
		return sim[:,end-hist_size+1:end]
	end
end

function semi_sync!(net::Net2,sim::sim_type,steps::Int,i_t0::Int,x...)
	nsigma_ins = rowvals(net.in_states)
	Nsize = length(net.nodes)
	buff = deepcopy(sim[:,i_t0])
	for t in i_t0:i_t0+steps-1
		n = 0
		for s in 1:Nsize
			s_indices = find(x->x.sync_order==s,net.nodes)
			for i in s_indices
				sim[i,t+1] = net.nodes[i].rule[buff[nsigma_ins[nzrange(net.in_states,i)]]]
				# sim[i,t+1] = net.nodes[i].rule[sim[nsigma_ins[nzrange(net.in_states,i)],t]]
				n += 1
			end
			buff[s_indices] = sim[s_indices,t+1]
			if n==Nsize
				break
			end
			# sim[:,t+1] = buff[:]
		end
	end
	return nothing
end

function fixed_delay!(net::Net2,sim::sim_type,steps::Int,i_t0::Int,x...)
	nsigma_ins = rowvals(net.in_states)
	# tau_ins = nonzeros(net.in_taus)
	# in_state_vals = nonzeros(net.in_states)
	# mu_in_state_vals = nonzeros(net.mu_in_states)
	for t in i_t0:i_t0+steps-1
		for (i,sigma) in enumerate(net.nodes)
			j_in = nzrange(net.in_states,i)
			input = [sim[nsigma_ins[nj],t-net.in_taus[nsigma_ins[nj],i]] for nj in j_in]
			sim[i,t+1] = sigma.rule[input]
			# mu_in_state_vals[j_in] = in_state_vals[j_in] = input
			# in_state_vals[j_in] = input
		end
	end
	return nothing
end

function smoothed_fixed_delay!(net::Net2,sim::sim_type,steps::Int,i_t0::Int,x...)
	nsigma_ins = rowvals(net.in_states)
	# tau_ins = nonzeros(net.in_taus)
	in_state_vals = nonzeros(net.in_states)
	mu_in_state_vals = nonzeros(net.mu_in_states)
	for t in i_t0:i_t0+steps-1
		for (i,sigma) in enumerate(net.nodes)
			j_in = nzrange(net.in_states,i)
			sim[i,t+1] = sigma.rule[in_state_vals[j_in]]
		end
		for (i,sigma) in enumerate(net.nodes)
			for j in nzrange(net.in_states,i)
				nsigma = nsigma_ins[j]
				tau = net.in_taus[nsigma,i]
				sigma_in = net.nodes[nsigma]
				if tau==0
					mu_in_state_vals[j] = in_state_vals[j] = sim[nsigma,t+1]
				else
					mu_in_state_vals[j] += (sim[nsigma,t+1]-sim[nsigma,t-tau])/(tau+1)
					in_state_vals[j] = sigma_in.state_range[findlast(mu_in_state_vals[j] .>= sigma_in.thresholds)+1]
				end
			end
		end
	end
	return nothing
end

function strob!(net::Net2,sim::sim_type,steps::Int,i_t0::Int,t_init::Int=0,x...)
	nsigma_ins = rowvals(net.in_states)
	tau_ins = nonzeros(net.in_taus)
	in_state_vals = nonzeros(net.in_states)
	for t in i_t0:i_t0+steps-1
		for (i,sigma) in enumerate(net.nodes)
			j_in = nzrange(net.in_states,i)
			for j in j_in
				if mod(t_init+t-i_t0+1,net.in_taus[nsigma_ins[j],i]+1)==0
					# sigma.sigma_in_state[nsigmak] = sim[sigmak,t-sigma.tau_in[nsigmak]]
					in_state_vals[j] = sim[nsigma_ins[j],t]
				end
			end
			sim[i,t+1] = sigma.rule[in_state_vals[j_in]]
		end
	end
	return nothing
end

function smoothed_strob!(net::Net2,sim::sim_type,steps::Int,i_t0::Int,t_init::Int=0,x...)
	nsigma_ins = rowvals(net.in_states)
	# tau_ins = nonzeros(net.in_taus)
	in_state_vals = nonzeros(net.in_states)
	mu_in_state_vals = nonzeros(net.mu_in_states)
	for t in i_t0:i_t0+steps-1
		for (i,sigma) in enumerate(net.nodes)
			j_in = nzrange(net.in_states,i)
			for j in j_in
				nsigma = nsigma_ins[j]
				tau = net.in_taus[nsigma,i]
				sigma_in = net.nodes[nsigma]
				if tau==0
					mu_in_state_vals[j] = in_state_vals[j] = sim[nsigma,t]
				elseif mod(t_init+t-i_t0+1,tau+1)==0
					mean_sigmak = mu_in_state_vals[j] = mean(sim[nsigma,t-tau:t])
					in_state_vals[j] = sigma_in.state_range[findlast(mean_sigmak.>= sigma_in.thresholds)+1]
				end
			end
			sim[i,t+1] = sigma.rule[in_state_vals[j_in]]
		end
	end
	return nothing
end

function obs_rand_cond!{T<:state_type}(net::Net2{T},steps::Int,obs_node::Union{AbstractString,Int,Vector{String},Vector{Int}},nrand_init::Int=prod(length(sigma.state_range) for sigma in net.nodes);out_dir::AbstractString="",file_tag::AbstractString="",update_mode::Symbol=:sy,noise_vector::Vector{Pair{Symbol,Float64}}=Pair{Symbol,Float64}[],init_h::Symbol=:default,forcing_rhythms::Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}=Dict{Symbol,Tuple{Vector{Bool},Vector{Int}}}(),init_constraints...)
	Nsize = length(net.nodes)
	dump_to_file = false
	if out_dir != ""
		if !isdir(out_dir)
			error(out_dir, " is not a valid directory to save results")
		elseif !isempty(readdir(out_dir))
			error(out_dir, " must be empty to save results")
		else
			dump_to_file = true
			println("Results will be saved into ",out_dir," folder")
		end
	end
	index_node_obs = Any[]
	if typeof(obs_node) <: Vector
		length(obs_node)!=0 || error("You must enter at least 1 node id or index")
		index_node_obs = typeof(obs_node)<:Vector{Int}? obs_node :  collect(findfirst(sigma->sigma.id==x,net.nodes) for x in obs_node)
		for (i,x) in enumerate(obs_node)
			0 < index_node_obs[i] <= Nsize || error("Unexisting node: $x")
		end
	else
		index_node_obs = typeof(obs_node)<:Int? obs_node : findfirst(sigma->sigma.id==obs_node,net.nodes)
		0 < index_node_obs <= Nsize || error("Unexisting node: $obs_node")
	end
	length_obs = length(index_node_obs)
	obs = length_obs >=2 ? Vector{Tuple{Vararg{T,length_obs}}}() : Vector{T}()
	state_space_size = prod(length(sigma.state_range) for sigma in net.nodes)
	if(state_space_size==0) && all(x->isinteger(x)&&x>0,length(sigma.state_range) for sigma in net.nodes)
		state_space_size = 10000000000
	end
	if nrand_init>state_space_size
        nrand_init = state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
	println("Calculating ",length_obs >=2 ? join((sigma.id for sigma in net.nodes[index_node_obs]),',') : net.nodes[index_node_obs].id ," dynamics over $nrand_init random initial conditions")
	i = 0
	hist_size = 1 + (update_mode == :sy ? 0 : maximum(net.in_taus))
	t_size = steps+hist_size
	sim_result = Array{T,2}(Nsize,t_size)
	sigma_ext_forcing = collect(sigma_symbol=>rhythm_values[2] for (sigma_symbol,rhythm_values) in forcing_rhythms)
	init_constraints = vcat(init_constraints,sigma_ext_forcing)
	if dump_to_file
		f = Vector{GZip.GZipStream}(length_obs)
		for (i,sigma_i) in enumerate(index_node_obs)
			filename = joinpath(out_dir,join([join([join([net.nodes[sigma_i].id,"T$t_size","r$nrand_init"],"_"),file_tag],""),"dat","gz"],"."))
			f[i] = GZip.open(filename,"a")
			# close(f[i])
		end
		println("Succesfully created files for storing simulations")
		# for (i,sigma_i) in enumerate(index_node_obs)
		# 	filename = joinpath(out_dir,join([join([join([net.nodes[sigma_i].id,"T$t_size","r$nrand_init"],"_"),file_tag],""),"dat","gz"],"."))
		# 	f[i] = GZip.open(filename,"a")
		# end
		if nrand_init<state_space_size
			while i < nrand_init
				rand_init_net!(net;update_mode=update_mode,init_h=init_h,init_constraints...)
				evolve_net!(net,steps;keep=false,update_mode=update_mode,sim=sim_result,noise_vector=noise_vector,forcing_rhythms=forcing_rhythms)
				for (fi,obs_i) in enumerate(index_node_obs)
					write(f[fi],map(Int8,sim_result[obs_i,:]'))
					# seek(f[fi],-1)
				end
				i += 1
				if mod(i,10000)==0
					@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
				end
			end
		else
			println("Exhaustive exploration of state space")
			node_indices = collect(vertices(net.graph))
			nodes_to_constraint = Int[]
			if !isempty(init_constraints)
				ids = [string(id) for (id,state) in init_constraints]
				allunique(ids) || error("Duplicate constraint specification")
				intersect(ids,(sigma.id for sigma in net.nodes))==ids || error("Trying to set unexisting node")
				all(i->i<:Union{Int,Bool,Vector{Int},Vector{Bool}},(typeof(state) for (id,state) in init_constraints)) || error("Wrong format to specify state")
				push!(nodes_to_constraint,(findfirst([sigma.id for sigma in net.nodes],id) for id in ids)...)
				all(map((x,y)->typeof(x)<:Vector?intersect(x,y)==x : intersect([x],y)==[x],(state for (id,state) in init_constraints),(sigma.state_range for sigma in net.nodes[nodes_to_constraint]))) || error("Trying to set node to a forbidden state")
				node_indices = setdiff(node_indices,nodes_to_constraint)
				states = Vector{Vector{Int}}(length(init_constraints))
				for (ni,nsigma) in enumerate(nodes_to_constraint)
					state = init_constraints[ni][2]
					states[ni] = isempty(state) ? [net.nodes[nsigma].def_state] : unique(state)
				end
			end
			ranges = Vector{Vector{Int}}(nv(net.graph))
			for nsigma in node_indices
				ranges[nsigma] = net.nodes[nsigma].state_range
			end
			for (ni,nsigma) in enumerate(nodes_to_constraint)
				ranges[nsigma] = states[ni]
			end
			itstspace = product(ranges...)
			for leave in itstspace
				init_net!(net;update_mode=update_mode,init_h=init_h,init_array=collect(leave))
				evolve_net!(net,steps;keep=false,update_mode=update_mode,sim_result=sim_result)
				for (fi,obs_i) in enumerate(index_node_obs)
					write(f[fi],map(Int8,sim_result[obs_i,:]'))
		#			seek(f[fi],-1)
				end
				#push!(obs,sim_result[index_node_obs,:]...)
				i += 1
				if mod(i,10000)==0
					@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
				end
			end
			nrand_init = i
		end
		status_files = zeros(Int,length_obs)
		for (fi,sigma_i) in enumerate(index_node_obs)
			status_files[fi] = close(f[fi])
		end
		if isempty(find(status_files))
			println("Simulations were succesfully stored for nodes ", join(obs_node,", "," and "))
		else
			println("Simulations couldn't be stored for nodes ", join(obs_node[find(status_files)],", "," nor "))
		end
		return status_files
	else
		if nrand_init<state_space_size
			while i < nrand_init
				rand_init_net!(net;update_mode=update_mode,init_h=init_h,init_constraints...)
				evolve_net!(net,steps;keep=false,update_mode=update_mode,sim=sim_result,noise_vector=noise_vector,forcing_rhythms=forcing_rhythms)
				if length_obs==1
					push!(obs,sim_result[index_node_obs,:]...)
				else
					push!(obs,(tuple(sim_result[index_node_obs,ti]...) for ti in 1:t_size)...)
				end
				i += 1
				if mod(i,10000)==0
					@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
				end
			end
		else
			println("Exhaustive exploration of state space")
			node_indices = collect(vertices(net.graph))
			nodes_to_constraint = Int[]
			if !isempty(init_constraints)
				ids = [string(id) for (id,state) in init_constraints]
				allunique(ids) || error("Duplicate constraint specification")
				intersect(ids,(sigma.id for sigma in net.nodes))==ids || error("Trying to set unexisting node")
				all(i->i<:Union{Int,Bool,Vector{Int},Vector{Bool}},(typeof(state) for (id,state) in init_constraints)) || error("Wrong format to specify state")
				push!(nodes_to_constraint,(findfirst([sigma.id for sigma in net.nodes],id) for id in ids)...)
				all(map((x,y)->typeof(x)<:Vector?intersect(x,y)==x : intersect([x],y)==[x],(state for (id,state) in init_constraints),(sigma.state_range for sigma in net.nodes[nodes_to_constraint]))) || error("Trying to set node to a forbidden state")
				node_indices = setdiff(node_indices,nodes_to_constraint)
				states = Vector{Vector{Int}}(length(init_constraints))
				for (ni,nsigma) in enumerate(nodes_to_constraint)
					state = init_constraints[ni][2]
					states[ni] = isempty(state) ? [net.nodes[nsigma].def_state] : unique(state)
				end
			end
			ranges = Vector{Vector{Int}}(nv(net.graph))
			for nsigma in node_indices
				ranges[nsigma] = net.nodes[nsigma].state_range
			end
			for (ni,nsigma) in enumerate(nodes_to_constraint)
				ranges[nsigma] = states[ni]
			end
			itstspace = product(ranges...)
			for leave in itstspace
				init_net!(net;update_mode=update_mode,init_h=init_h,init_array=collect(leave))
				evolve_net!(net,steps;keep=false,update_mode=update_mode,sim_result=sim_result)
				if length_obs==1
					push!(obs,sim_result[index_node_obs,:]...)
				else
					push!(obs,(tuple(sim_result[index_node_obs,ti]...) for ti in 1:t_size)...)
				end
				#push!(obs,sim_result[index_node_obs,:]...)
				i += 1
				if mod(i,10000)==0
					@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
				end
			end
			nrand_init = i
		end
		return reshape(obs,round(Int,length(obs)/nrand_init),nrand_init)
	end
end
