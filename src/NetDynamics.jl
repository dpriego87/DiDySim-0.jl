using StatsBase

function init_node!(sigma::BNode,state::Union{Tuple{Vararg{Integer}},Integer};init_history::Symbol=:default)
	if typeof(state)<:Tuple
		for state_i in state
			!isempty(filter(x->x==state_i,sigma.state_range)) || error("Trying to set node $(sigma.id) to a forbidden state, value $state_i is out of possible states range $(sigma.state_range)")
		end
	else
		!isempty(filter(x->x==state,sigma.state_range)) || error("Trying to set node $(sigma.id) to a forbidden state, value $state is out of possible states range $(sigma.state_range)")
	end
	if init_history==:default
		fill!(sigma.state,sigma.def_state)
	elseif init_history==:rand
		sigma.state[:]=rand(sigma.state_range,length(sigma.state))
	elseif init_history!=:last
		error("Invalid initialization flag. Available options are\n:default for default, :rand for random")
	end
	sigma.state[end-length(state)+1:end] = collect(state)
	return nothing
end

function init_net!(netstate::Vector{BNode};init_t0::Symbol=:custom,init_history::Symbol=:last,init_vector::Vector{Int}=Int[],init_constraints...)
	if init_t0==:default
		for sigma in netstate
			init_node!(sigma,sigma.def_state;init_history=init_history)
		end
	elseif init_t0==:rand
		for sigma in netstate
			init_node!(sigma,rand(sigma.state_range);init_history=init_history)
		end
	elseif init_t0==:custom
		isempty(init_vector)||length(init_vector)==length(netstate)||error("Initial configuration vector must have the same length of network size")
		for (n,sigma) in enumerate(netstate)
			init_node!(sigma,init_vector[n];init_history=init_history)
		end
	else
		error("Invalid initialization flag. Available options are\n:default for default, :rand for random :custom start from last available configuration")
	end

	for (id,state) in init_constraints
		ni = findfirst(collect(sigma.id for sigma in netstate),string(id))
		ni!=0 || error("Trying to set unexisting $id node")
		init_node!(netstate[ni],state,init_history=:last)
	end
	for sigma in netstate
		sigma.sigma_in_state[:] = [netstate[sigmak].state[end-sigma.tau_in[nsigmak]] for (nsigmak,sigmak) in enumerate(sigma.sigma_in)]
	end
	for sigma in netstate
		for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
			if sigma.tau_in[nsigmak]==0
				sigma.sigma_in_state[nsigmak] = sigma.sigma_in_state_mu[nsigmak] = netstate[sigmak].state[end]
			else
				mean_sigmak = sigma.sigma_in_state_mu[nsigmak] = mean(netstate[sigmak].state[end-sigma.tau_in[nsigmak] : end])
				sigma.sigma_in_state[nsigmak] = netstate[sigmak].state_range[findlast(mean_sigmak .>= netstate[sigmak].thresholds)+1]
			end
		end
	end
	return hcat((sigma.state for sigma in netstate)...)'
end

function evolve_net!(netstate::Vector{BNode},steps::Integer;keep::Bool=true,update_mode::Symbol=:sy,sim::Array{Int8,2}=Array{Int8}(length(netstate),steps + 1 + (update_mode == :sy ? 0 : maximum(vcat([sigma.tau_in for sigma in netstate]...)))),t0::Integer=0)
	update_schemes = Dict(:sy=>sync!,:fd=>fixed_delay!,:md=>smoothed_fixed_delay!,:st=>strob!,:mst=>smoothed_strob!,:asyr=>async_rand!,:asyrp=>async_rand_perm!)
	haskey(update_schemes,update_mode) || error("Unexisting update scheme. Available schemes are :sy,:fd,:md,:sd,:smd")
	steps>=0 || error("Steps must be a positive integer")
	if update_mode == :sy
		tau_max = 0
	else
		tau_max = maximum(vcat([sigma.tau_in for sigma in netstate]...))
	end
	i_t0 = tau_max + 1
	# sim = Array{Int8}(length(netstate),steps + i_t0)
	# Initialize simulation table (from t0 - tau0 to t0)
	size(sim) == (length(netstate),steps + i_t0) || error("Array for storing results has incorrect dimensions $(size(sim)). Expected array size is ($(length(netstate)),$(steps + i_t0))")
	if update_mode == :sy
		sim[:,1] = [sigma.state[end] for sigma in netstate]
	else
		# for (nsigma,sigma) in enumerate(netstate)
		# 	tau0 = i_t0 - length(sigma.state) + 1
		# 	sim[nsigma,tau0:i_t0] = sigma.state
		# 	sim[nsigma,1:tau0-1] = fill(sigma.state[1],tau0-1)
		# end
		for (nsigma,sigma) in enumerate(netstate)
			sim[nsigma,1:i_t0] = sigma.state
		end
	end
	if update_mode == :st || update_mode == :mst
		update_schemes[update_mode](netstate,sim,steps,i_t0,t0)
	else
		update_schemes[update_mode](netstate,sim,steps,i_t0)
	end
	for (nsigma,sigma) in enumerate(netstate)
		sigma.state[end-tau_max:end] = sim[nsigma,end-tau_max:end]
	end
	if keep
		result = sim
	else
		result = sim[:,end-i_t0+1:end]
	end
	return result
end

function sync!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer)
	for t in 1:steps
		sim[:,t+1] = [sigma.F[tuple([sim[sigmak,t] for sigmak in sigma.sigma_in]...)] for sigma in netstate]
	end
	nothing
end

function async_rand!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer)
	for t in 1:steps
		sim[:,t+1] = sim[:,t]
		i = rand(eachindex(netstate))
		sim[i,t+1] = netstate[i].F[tuple([sim[sigmak,t] for sigmak in netstate[i].sigma_in]...)]
	end
	nothing
end

function async_rand_perm!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer)
	for t in 1:steps
		for i in shuffle(eachindex(netstate))
			sim[i,t+1] = netstate[i].F[tuple([sim[sigmak,t] for sigmak in netstate[i].sigma_in]...)]
		end
	end
	nothing
end

function fixed_delay!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer)
	for t in i_t0:i_t0+steps-1
		sim[:,t+1] = [sigma.F[tuple([sim[sigmak,t-sigma.tau_in[nsigmak]] for (nsigmak,sigmak) in enumerate(sigma.sigma_in)]...)] for sigma in netstate]
	end
	nothing
end

function smoothed_fixed_delay!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer)
	for (nsigma,sigma) in enumerate(netstate)
		for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
			if sigma.tau_in[nsigmak]==0
				sigma.sigma_in_state_mu[nsigmak] = sigma.sigma_in_state[nsigmak] = sim[sigmak,i_t0]
 			else
				mean_sigmak = mean(sim[sigmak,i_t0-sigma.tau_in[nsigmak] : i_t0])
				sigma.sigma_in_state_mu[nsigmak] = mean_sigmak
				sigma.sigma_in_state[nsigmak] = netstate[sigmak].state_range[findlast(mean_sigmak .>= netstate[sigmak].thresholds)+1]
			end
		end
	end
	for t in i_t0:i_t0+steps-1
		sim[:,t+1] = [sigma.F[tuple(sigma.sigma_in_state...)] for sigma in netstate]
		for (nsigma,sigma) in enumerate(netstate)
			for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
				if sigma.tau_in[nsigmak]==0
					sigma.sigma_in_state_mu[nsigmak] = sigma.sigma_in_state[nsigmak] = sim[sigmak,t+1]
				else
					sigma.sigma_in_state_mu[nsigmak] += (sim[sigmak,t+1]-sim[sigmak,t-sigma.tau_in[nsigmak]])/(sigma.tau_in[nsigmak]+1)
					sigma.sigma_in_state[nsigmak] = netstate[sigmak].state_range[findlast(sigma.sigma_in_state_mu[nsigmak] .>= netstate[sigmak].thresholds)+1]
				end
			end
		end
	end
	nothing
end

function strob!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer,t0::Integer=0)
	for t in i_t0:i_t0+steps-1
		for (nsigma,sigma) in enumerate(netstate)
			for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
				if mod(t0+t-i_t0+1,sigma.tau_in[nsigmak]+1)==0
					# sigma.sigma_in_state[nsigmak] = sim[sigmak,t-sigma.tau_in[nsigmak]]
					sigma.sigma_in_state[nsigmak] = sim[sigmak,t]
				end
			end
			sim[nsigma,t+1] = sigma.F[tuple(sigma.sigma_in_state...)]
		end
	end
	nothing
end

function smoothed_strob!(netstate::Vector{BNode},sim::Array{Int8,2},steps::Integer,i_t0::Integer,t0::Integer=0)
	# for (nsigma,sigma) in enumerate(netstate)
	# 	for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
	# 		if sigma.tau_in[nsigmak]==0
	# 			sigma.sigma_in_state[nsigmak] = sigma.sigma_in_state_mu[nsigmak] = sim[sigmak,i_t0]
	# 		else
	# 			mean_sigmak = sigma.sigma_in_state_mu[nsigmak] = mean(sim[sigmak,i_t0-sigma.tau_in[nsigmak] : i_t0])
	# 			sigma.sigma_in_state[nsigmak] = netstate[sigmak].state_range[findlast(mean_sigmak .>= netstate[sigmak].thresholds)+1]
	# 		end
	# 	end
	# end
	for t in i_t0:i_t0+steps-1
		for (nsigma,sigma) in enumerate(netstate)
			for (nsigmak,sigmak) in enumerate(sigma.sigma_in)
				if sigma.tau_in[nsigmak]==0
					sigma.sigma_in_state[nsigmak] = sigma.sigma_in_state_mu[nsigmak] = sim[sigmak,t]
				elseif mod(t0+t-i_t0+1,sigma.tau_in[nsigmak]+1)==0
					mean_sigmak = sigma.sigma_in_state_mu[nsigmak] = mean(sim[sigmak,t-sigma.tau_in[nsigmak] : t])
					sigma.sigma_in_state[nsigmak] = netstate[sigmak].state_range[findlast(mean_sigmak .>= netstate[sigmak].thresholds)+1]
				end
			end
			sim[nsigma,t+1] = sigma.F[tuple(sigma.sigma_in_state...)]
		end
	end
	nothing
end

function obs_rand_cond!(netstate::Vector{BNode},steps::Integer,obs_node_id::AbstractString,nrand_init::Integer=prod([length(sigma.state_range) for sigma in netstate]);update_mode::Symbol=:sy,init_constraints...)
	index_node_obs = findfirst(sigma->sigma.id==obs_node_id,netstate)
	state_space_size = prod([length(sigma.state_range) for sigma in netstate])
    if nrand_init>state_space_size
        nrand_init = state_space_size
        println("Oversampling error: Random initial conditions bigger than possible number of configurations")
    end
	if index_node_obs==0
		error("Unexisting node: $obs_node_id")
	end
	obs = Int8[]
	println("Calculating averaged $(netstate[index_node_obs].id) simulation over $nrand_init random initial conditions")
	tags = Dict{Tuple{Vararg{Int8}},Bool}()
	i = 0
	tau_max = maximum([maximum(sigma.tau_in) for sigma in netstate])
	sim_result = update_mode == :sy ? Array{Int8}(length(netstate),steps+1) : Array{Int8}(length(netstate),steps+tau_max+1)
	if nrand_init<state_space_size
		while i < nrand_init
			ri = init_net!(netstate;init_t0=:rand,init_constraints...)
			if haskey(tags,tuple(ri...))
				#println("Configuration already visitted")
				continue
			end
			tags[ri...] = true
			evolve_net!(netstate,steps,update_mode=update_mode,sim=sim_result)
			push!(obs,sim_result[index_node_obs,:]...)
			i += 1
			if mod(i,10000)==0
				@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
			end
		end
	else
		net_size = length(netstate)
		state_ranges = [sigma.state_range for sigma in netstate]
		itstspace = product(state_ranges...)
		index_to_constraint = Dict{Integer,Integer}()
		for (id,state) in init_constraints
			node_to_init = findfirst(sigma->sigma.id==string(id),netstate)
			node_to_init!=0 || error("Trying to set unexisting $id node")
			index_to_constraint[node_to_init] = state
		end
		for leave in itstspace
			jump = false
	        for k in eachindex(netstate)
				init_node!(netstate[k],leave[k])
	        end
			for (i_const,state) in index_to_constraint
				if leave[i_const]!=state
					jump = true
					break
				end
			end
			if jump
				continue
			end
			evolve_net!(netstate,steps,update_mode=update_mode,sim=sim_result)
			push!(obs,sim_result[index_node_obs,:]...)
			i += 1
			if mod(i,10000)==0
				@printf("%0.2f%% progress, %d out of %d initial conditions explored\n", 100*i/nrand_init, i, nrand_init)
			end
		end
		nrand_init = i
	end
	return reshape(obs,round(Int,length(obs)/nrand_init),nrand_init)
end
