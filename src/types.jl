reg_function = Union{Function,Dict{Vector{Int8},Int8},Dict{BitArray{1},Bool}}
state_type = Union{Bool,Int8}
sim_type = Union{BitArray{2},Array{Int8,2},Array{Bool,2}}
mutable struct BNode
    id::AbstractString               # Name
    state::Vector{Int8}              # Current state
    def_state::Int8              # Default state
    state_range::Vector{Int8}    # Possible states
    sigma_in::Vector{Int16}
    sigma_in_state::Vector{Int8}
    sigma_in_state_mu::Vector{Float16}
    sigma_out::Vector{Int}
    F::Dict{Tuple{Vararg{Int8}},Int8}
    tau_in::Vector{UInt16}
    tau_out::UInt16
    thresholds::Vector{Float16}
    function BNode()
        this = new()
        this.state = [false]
        this.def_state = false
        this.state_range = [false,true]
        this.sigma_in = Vector{Int16}()
        this.sigma_in_state = Vector{Int8}()
        this.sigma_in_state_mu = Vector{Float16}()
        this.sigma_out = Vector{Int}()
        this.F = Dict{Tuple{Vararg{Int8}},Int8}()
        this.tau_in = Vector{UInt16}()
        this.tau_out = 1
        this.thresholds = Vector{Float16}()
        return this
    end
    function BNode(id::AbstractString)
        this = new()
        this.id = id
        this.state = [false]
        this.def_state = false
        this.state_range = [false,true]
        this.sigma_in = Vector{Int16}()
        this.sigma_in_state = Vector{Int8}()
        this.sigma_in_state_mu = Vector{Float16}()
        this.sigma_out = Vector{Int}()
        this.F = Dict{Tuple{Vararg{Int8}},Int8}()
        this.tau_in = Vector{UInt16}()
        this.tau_out = 1
        this.thresholds = Vector{Float16}()
        return this
    end
    #Node(id::AbstractString,state::Integer,def_state::Integer,state_range::Vector{Integer},sigma_in::Vector{Node})=new(id,state,def_state,state_range,sigma_in,Dict{Tuple{Vararg{Int8}},Integer}())
    #Node(id::AbstractString,state::Integer,def_state::Integer,state_range::Vector{Integer},sigma_in::Vector{Node},F::Dict{Tuple{Vararg{Int8}},Integer})=new(id,state,def_state,state_range,sigma_in,F)
end


# type BNode{T<:Union{Int8,AbstractString}}
#     id::T                           # Name
#     state::Int8                  # Current state
#     def_state::Int8              # Default state
#     state_range::Vector{Int8}    # Possible states
#     sigma_in::Vector{BNode{T}}
#     F::Dict{Tuple{Vararg{Int8}},Int8}
#     BNode() = new(false,false,[false,true],Vector{BNode{T}}(),Dict{Tuple{Vararg{Integer}},Integer}())
# end

mutable struct Net
    state::Vector{BNode}
    graph::DiGraph
    Net(state::Vector{BNode},graph::DiGraph) = new(state,graph)
end

mutable struct BNode2{T<:state_type}
    id::AbstractString               # Name
    def_state::T             # Default state
    state_range::Vector{T}    # Possible states
    rule::reg_function
    thresholds::Vector{Float16}
    sync_order::Int
    function BNode2{T}(id::AbstractString,def_state::Union{Int,Bool},state_range::Union{Vector{Int},Vector{Bool}}) where {T}
        vectorType = T==Bool ? BitArray{1} : Vector{T}
        this = new()
        this.id = id
        this.def_state = def_state
        this.state_range = state_range
        this.rule = Dict{vectorType,T}()::reg_function
        this.thresholds = Vector{Float16}()
        this.sync_order = 1
        return this
    end
    function BNode2{T}(id::AbstractString,def_state::Union{Int,Bool},state_range::Union{Vector{Int},Vector{Bool}},sync_order::Int) where {T}
        vectorType = T==Bool ? BitArray{1} : Vector{T}
        this = new()
        this.id = id
        this.def_state = def_state
        this.state_range = state_range
        this.rule = Dict{vectorType,T}()::reg_function
        this.thresholds = Vector{Float16}()
        this.sync_order = sync_order
        return this
    end
end

mutable struct Net2{T<:state_type}
    nodes::Vector{BNode2{T}}
    graph::DiGraph
    in_taus::SparseMatrixCSC{UInt16,Int}
    in_states::SparseMatrixCSC{T,Int}
    mu_in_states::SparseMatrixCSC{Float16,Int}
    sim_buffer::Union{BitArray{2},Array{T,2}}
    Net2{T}(nodes::Vector{BNode2{T}},graph::DiGraph) where {T<:state_type} = new(nodes,graph,spzeros(nv(graph),nv(graph)),adjacency_matrix(graph),adjacency_matrix(graph),T==Bool ? falses(nv(graph),1) : zeros(T,nv(graph),1))
    Net2{T}(nodes::Vector{BNode2{T}},graph::DiGraph,in_taus::SparseMatrixCSC{Int,Int}) where {T<:state_type}= (nv(graph),nv(graph)) == size(in_taus) ? new(nodes,graph,in_taus,adjacency_matrix(graph),adjacency_matrix(graph),T==Bool ? falses(nv(graph),maximum(in_taus)+1) : zeros(T,nv(graph),maximum(in_taus)+1)) : error("tau_matrix must have the same dimensions as adjacency matrix of network")
end

    #Node(id::AbstractString,state::Integer,def_state::Integer,state_range::Vector{Integer},sigma_in::Vector{Node})=new(id,state,def_state,state_range,sigma_in,Dict{Tuple{Vararg{Int8}},Integer}())
    #Node(id::AbstractString,state::Integer,def_state::Integer,state_range::Vector{Integer},sigma_in::Vector{Node},F::Dict{Tuple{Vararg{Int8}},Integer})=new(id,state,def_state,state_range,sigma_in,F)
