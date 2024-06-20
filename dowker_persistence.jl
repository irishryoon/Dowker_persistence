"""
Computes Dowker persistence with F2 coefficients. 

@author Iris Yoon
irishryoon@gmail.com
"""
module Dowker

include("Eirene_var.jl")
using Combinatorics
using Distances
using .Eirene_var
using Interact
using IJulia
using JLD
using Measures
using Plots
using Plots.PlotMeasures
using Random
using Statistics
using StatsBase
using Printf
export 
    compute_Witness_persistence,
    find_Dowker_cycle_correspondence,
    compute_distance,
    plot_barcode,
    select_intervals,
    get_Witness_cyclerep,
    plot_PD,
    plot_P_Q,
    plot_cycle,
    get_1simplices,
    get_2simplices,
    plot_Dowker_complex
plotly()

#############################################################################################
# MAIN FUNCTION 
#############################################################################################


#################################################################################
# functions for visualizations
#################################################################################

function plot_barcode(barcode; 
    color = :grey56, # default bar color
    selected_bars = [], # index of bars to highlight
    epsilon = missing, # if provided, only highlight the portion of bars on the right side of epsilon
    selection_color = :deeppink2,  # highlight color
    v_line = [], # if provided, draw vertical lines at values
    return_perm = false, # whether to return the permutation index or not
    kwargs...)

    # adjust linewidth, depending on number of intervals
    n = size(barcode,1)

    # find ordering according to birth time
    perm = sortperm(barcode[:,1])

    # non-inf maximum death time
    if filter(!isinf,barcode[:,2]) != []
        death_max = maximum(filter(!isinf,barcode[:,2])) 
    else
        death_max = maximum(barcode[:,1]) * 2
    end

    p = plot(framestyle = :box,
            top_margin = 5 * Plots.PlotMeasures.mm, 
            bottom_margin = 5 * Plots.PlotMeasures.mm, 
            yaxis = nothing;
            kwargs...)
    
    # plot all bars
    idx = 1
    for i in perm
        birth = barcode[i,1]
        death = barcode[i,2]
        
        # assign a death time if bar has infinite death time 
        if isinf(death)
            death = death_max * 1.2
        end
        if i in selected_bars
            
            # if epsilon is missing, highlight the entire bar
            if ismissing(epsilon)
                plot!(p,[birth, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
            
            # if epsilon is provided, only highlight the portion of the bar on the right side of epsilon    
            else 
                if birth <= epsilon
                    plot!(p,[birth, epsilon], [idx, idx], legend = false, linecolor = color, hover = "class " *string(i); kwargs...)
                    plot!(p,[epsilon, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
                else
                    plot!(p,[birth, death], [idx, idx], legend = false, linecolor = selection_color, hover = "class " *string(i); kwargs...)
                end
            end
        else
            plot!(p,[birth,death],[idx,idx], legend = false, linecolor = color, hover = "class " *string(i); kwargs...)
        end
        idx += 1
    end

    # plot vertical lines 
    if v_line != []
        plot!(v_line, seriestype="vline", linestyle = :dot, linecolor = :red)
    end

    ylims!((-1, n+1))
    
    if return_perm == true
        return p, perm
    else
        return p
    end
end



"""
    plot_cycle_single()
Plots a single point cloud P in 2-dimensions and a 1-dimensional cycle. 
"""
function plot_cycle_single(P; cycle = [], cycle_color = "black", cycle_lw = 5, kwargs...)
    # P: array of size (2, n) or (3,n)
    # cycle: [[v1, v2], [v3, v4], ...  ]
    
    # plot points P
    p = plot(P[1,:], P[2,:], 
            seriestype = :scatter, 
            label = "",
            framestyle = :box,
            xaxis = nothing,
            yaxis = nothing;
            kwargs...)
    
    # plot cycle
    for item in cycle
        p1, p2 = item
        p1_x, p1_y = P[:,p1]
        p2_x, p2_y = P[:,p2]
        plot!(p, [p1_x, p2_x], [p1_y, p2_y], color = cycle_color, lw = cycle_lw, label ="")
    end
    
    return p
end



"""
    plot_P_Q()
Plots point clouds P and Q in 2-dimensions
"""
function plot_P_Q(P, # array of size (m, 2)
                  Q; # array of size (n, 2)
                  P_color = "#008181", 
                  P_label = "P",
                  P_markersize = 5,
                  P_marker = :circle,
                  Q_color = "#ff8d00",
                  Q_label = "Q",
                  Q_markersize = 5,
                  Q_marker = :xcross,
                  kwargs...)
   # plot points P and Q on a square
    
    p = plot(framestyle = :box, yaxis = nothing, xaxis = nothing; kwargs...)
    
    # plot P
    scatter!(p, P[:,1], P[:,2], color = P_color, label = P_label, markersize = P_markersize, marker = P_marker)
    
    # plot Q
    scatter!(p, Q[:,1], Q[:,2], color = Q_color, label = Q_label, markersize = Q_markersize, marker = Q_marker)
    return p
end


"""
    plot_cycle
Plots both points P, Q and a 1-dimensional cycle. User can specify whetheer the cycle exists among P or Q.
"""
function plot_cycle(P, # array of size (m,2)
                    Q; # array of size (n,2)
                  cycle = [],
                  cycle_loc = "P",
                  cycle_color = :deeppink,
                  cycle_linewidth = 5,
                  P_color = "#008181", 
                  P_label = "P",
                  P_markersize = 5,
                  P_marker = :circle,
                  Q_color = "#ff8d00",
                  Q_label = "Q",
                  Q_markersize = 5,
                  Q_marker = :xcross,
                  kwargs...)
# plot one-dimensional cycle 
        
    p = plot_P_Q(P, Q, 
                 P_color = P_color, P_label = P_label, P_markersize = P_markersize, P_marker = P_marker,
                 Q_color = Q_color, Q_label = Q_label, Q_markersize = Q_markersize, Q_marker = Q_marker;
                 kwargs...)

    # specifiy the 0-simplices
    if cycle_loc == "P"
        PC = P
    else
        PC = Q
    end
    
    for simplex in cycle
        v1, v2 = simplex 
        v1_theta, v1_phi = PC[v1,1], PC[v1,2]
        v2_theta, v2_phi = PC[v2,1], PC[v2,2]
        
        plot!(p, [v1_theta, v2_theta], [v1_phi, v2_phi], label = "", color = cycle_color, lw = cycle_linewidth)

    end
    return p
end


#################################################################################
# Functions implementing Dowker's Theorem
#################################################################################
"""
    apply_Dowker(W_PQ, W_QP; <keyword arguments>)
Given Witness filtrations `W_PQ` (landmark: P, witness: Q) and `W_QP` (landmark: Q, witness: P), 
use Dowker's Theorem to find the correspondence between bars in barcode(W(P,Q)) and barcode(W(Q,P)).

### Arguments
- `W_PQ`(dict): Output of `compute_Witness_persistence(D_P_Q, maxdim = dim)`
- `W_QP`(dict): Output of `compute_Witness_persistence(D_Q_P, maxdim = dim)`
- dim(int): dimension. Defaults to 1

### Outputs
- `P_to_Q` (dict): of correspondence between barcode(W(P,Q)) and barcode(W(Q,P)).
        `P_to_Q[i] = j` implies that bar `i` in barcode(W(P,Q)) matches bar `j` in barcode(W(Q,P))
"""
function apply_Dowker(
    W_PQ,
    W_QP;
    dim = 1)

    P_to_Q = Dict()
    barcode_W_PQ = barcode(W_PQ["eirene_output"], dim = dim)
    barcode_W_QP = barcode(W_QP["eirene_output"], dim = dim)
    n = size(barcode_W_PQ)[1]
    for i=1:n
        birth_rows = findall(x->x==1, (barcode_W_QP[:, 1] .== barcode_W_PQ[i,1]))
        if size(birth_rows)[1] > 1
            print("ERROR: multiple bars with same birth time")
            break
        else
            P_to_Q[i] = birth_rows[1]
        end
    end
    return P_to_Q
end


"""
    find_barycentric_subdivision(simplex)
Finds the barycentric subdivision of a 1-dimensional simplex

### Arguments
- simplex: 1-dimensional simplex of form [i, j], consisting of the i-th and j-th vertices

### Returns
- a list containing the vertices of the barycentric subdivision
"""
function find_barycentric_subdivision(simplex)
    return [[simplex[1]], simplex, [simplex[2]]]
    
end

"""
    find_witness_column(relations, rows)
Given rows in a relations matrix, return a witness column. In particular, return witness column with smallest index

### Arguments
- relations: (array) binary relations matrix
- rows: (list) of rows

### Returns
- col: (int) that witnesses the given rows
"""
function find_witness_column(relations, rows)
    sub_relations = relations[rows, :]
    n_cols = size(relations, 2)
    n_rows = size(rows, 1)

    # compute the sum of all rows in sub matrix
    S = sum(sub_relations, dims = 1)
    idx = findfirst(x -> x == n_rows, S)
    
    if idx == nothing
        print("There is no witness for selected rows")
    else
        col = idx[2]
        return col
    end
end

"""
    barycentric_to_column_complex(barycentric_simplex, relations)
Given a list of vertices in the barycentric subdivision, return 1-simplices in the column complex it maps to.
That is, given a simplex in the barycentric subdivision, return all witnesses of the vertices.

### Arguments
- relations: (array) relations matrix

### Returns
"""
function barycentric_to_column_complex(barycentric_simplex, relations)
    
    mapped_cols = [find_witness_column(relations, vertex) for vertex in barycentric_simplex]
    mapped_1simplices = [[mapped_cols[1], mapped_cols[2]], [mapped_cols[2], mapped_cols[3]]]
    mapped_1simplices = [item for item in mapped_1simplices if item[1] != item[2]]
    
    return mapped_1simplices
end

"""
    find_column_cycle_via_Dowker(row_cycle, relations)
Given a cycle in the row complex, find the corresponding cycle in the column complex using Dowker's Theorem
"""
function find_column_cycle_via_Dowker(row_cycle, relations)
    # check that input is a cycle
    col_cycle = []
    for row_simplex in row_cycle
        barycentric = find_barycentric_subdivision(row_simplex)
        col_simplex = barycentric_to_column_complex(barycentric, relations)
        append!(col_cycle, col_simplex)
    end
    col_cycle = [sort(item) for item in col_cycle]
    return col_cycle
end



"""
    find_Dowker_cycle_correspondence(cycle_W_PQ, param, D_P_Q)
Given a cycle in the Witness complex W(P,Q) and a parameter at which the cycle is nontrivial, 
find its corresponding cycle in W(Q,P) via Dowker's Theorem. 
"""
function find_Dowker_cycle_correspondence(cycle_W_PQ, param, D_P_Q)
    # find the binary relations matrix at parameter epsilon
    relations = replace(x -> x <= param, D_P_Q)

    # find corresponding cycle in W(Q,P) using Dowker Theorem
    cycle_W_QP = find_column_cycle_via_Dowker(cycle_W_PQ, relations)
    
    return cycle_W_QP
end



#################################################################################
# Functions related to Witness complexes
#################################################################################

"""
    select_vertices_of_Witness_complex(D, threshold)
Given a cross-distance matrix D, select the vertices that exist at W^{threshold}. Since W^{threshold} may have fewer vertices than given (rows of D), we need to index the vertices of W^{threshold}.

### Arguments
- `D`: (array) of dross-distance matrix
- `threshold`: (float) parameter to build W^{threshold}

### Returns
- if all vertices are selected, return nothing, nothing
- if a subset of vertices are selected, then return the following
    - `Wpsi_vertex_to_default_vertex`: (dict)
        `Wpsi_vertex_to_default_vertex[i] = j` means vertex i of Wpsi corresponds to vertex j in default vertex (or row j in matrix D)
    - `default_vertex_to_Wpsi_vertex`: (dict) 
        `default_vertex_to_Wpsi_vertex[i] = j` means vertex i (or row i of matrix D) corresponds to vertex j in Wpsi
"""
function select_vertices_of_Witness_complex(D; threshold = Inf)
    
    # Find rows of D that has at least one entry <= threshold
    D_binary = replace(x -> x <= threshold, D)
    idx = findall(x -> x >0 , [(sum(D_binary, dims = 2)...)...])
    

    # Need to index the vertices of Wpsi
   
    if size(idx, 1) == size(D, 1)
        return nothing, nothing
    else
        # create correspondence between vertices of Wepsilon and default vertices corresponding to rows of D
        Wpsi_vertex_to_default_vertex = Dict(i => idx[i] for i=1:size(idx,1))
        default_vertex_to_Wpsi_vertex = Dict(val => key for (key, val) in Wpsi_vertex_to_default_vertex)
    
        return Wpsi_vertex_to_default_vertex, default_vertex_to_Wpsi_vertex
    end
end


function build_Witness_complex(
    D::Array, 
    threshold::Float64;
    maxdim::Int64=2)
    
    # select vertices that exist at given threshold
    Wepsilon_vertex_to_default_vertex, default_vertex_to_Wepsilon_vertex = select_vertices_of_Witness_complex(D, threshold = threshold)    
    if Wepsilon_vertex_to_default_vertex == nothing
        D_sub = D
        n = size(D, 1)
    else
        W_idx = sort(collect(keys(default_vertex_to_Wepsilon_vertex)))
        D_sub = D[W_idx, :]
        n = size(W_idx,1)
    end
        
    Witness_complex = []
    fv = []
    for d = 0:maxdim
        simplices = []
        fv_dim = []

        candidates = collect(combinations(1:n, d+1))
        for item in candidates
            t = minimum(maximum(D_sub[item,:], dims = 1))
            if t <= threshold
                push!(simplices, item)
                push!(fv_dim, t)
            end
        end

        push!(Witness_complex, simplices)
        push!(fv, fv_dim)
    end
    
    # turn the list of simplices to dictionary
    # Create dictionaries
    W_simplex2index = Dict([i] => i for i =1:n)

    # higher dimensional simplices
    for d = 1:maxdim
        for (idx, val) in enumerate(Witness_complex[d+1])
            W_simplex2index[val] = idx
        end
    end
    
    # rverse simplex2index
    W_index2simplex = Dict((val, size(key,1)-1) => key for (key, val) in W_simplex2index)
    
    return W_index2simplex, W_simplex2index, fv, Wepsilon_vertex_to_default_vertex, default_vertex_to_Wepsilon_vertex
end


"""
    compute_Witness_persistence(D; <keyword arguments>)
Given a cross-distance matrix `D`, build the Witness filtration using rows as landmarks (0-simplices) and columns as witnesses.

### Arguments
- `D`: Distance matrix between landmarks and witnesses.
    rows:landmarks
    columns: witnesses
- `maxdim`: Maximum dimension of interest. Defaults to 1
- `param_max`: maximum parameter to build the Witness filtration.

### Outputs 
- `W`: Dictionary containing the following information.
    - `param_max`: threshold value used in `build_Witness_complex`
    - `index2simplex`: indexing of simplices in the Witness complex
    - `simpelx2index`: referse of index2simplex
    - `distance_matrix`: D 
    - `fv`: birth times of simplices in the Witness complex
    - `eirene_output`: dictionary output of Eirene on the Witness filtration. 
    - `W_vertex_to_default_vertex`: Output of `select_vertices_of_Witness_complex``
    - `default_vertex_to_W_vertex`: Output of `select_vertices_of_Witness_complex`
"""
function compute_Witness_persistence(
    D::Array;
    maxdim::Int64 = 1,
    param_max = Inf)
    
    # select max parameter
    if param_max == Inf
        param_max = minimum(maximum(D, dims = 1))
    end
        
    # Index simplices in the Witness complex and find face value (birth time) in filtration.
    W_index2simplex, W_simplex2index, W_fv, W_vertex_to_default_vertex, default_vertex_to_W_vertex = build_Witness_complex(D, param_max, maxdim = maxdim + 1)
    
    # prepare input for Eirene
    W_rv, W_cp, W_ev, W_fv = create_CSC_input_from_fv(W_fv, W_index2simplex, W_simplex2index)
    
    # run Eirene
    C_W = eirene(rv = W_rv, cp = W_cp, ev = W_ev, fv = W_fv, record = "all", maxdim = maxdim)
    
    # save relevant info
    W = Dict()
    W["param_max"] = param_max
    W["index2simplex"] = W_index2simplex
    W["simplex2index"] = W_simplex2index
    W["distance_matrix"] = D
    W["fv"] = W_fv
    W["eirene_output"] = C_W
    W["W_vertex_to_default_vertex"] = W_vertex_to_default_vertex
    W["default_vertex_to_W_vertex"] = default_vertex_to_W_vertex
    
    return W
end

"""
    get_Witness_cyclerep(W, class_num)
Find the classrep of a Witness filtration in vertex format.
"""
function get_Witness_cyclerep(W;class_num = nothing, dim = 1)
    
    cyclerep_idx = classrep(W["eirene_output"], dim = dim, class = class_num)
    
    cyclerep = []
    for i in cyclerep_idx
       push!(cyclerep, W["index2simplex"][(i, dim)]) 
    end
    return cyclerep
end
        

########################################################################################
# functions for working with Eirene's CSC format
########################################################################################
function create_CP(
    ev)
    # create `cp` for CSC format input to Eirene 
    """
    --- input ---
    ev: (array)
        ev[d]: number of (d-1) dimensional simplices
    --- output ---
    cp: (array) colptr to be used as input to Eirene
    """
    cp = []
    for d = 1:size(ev,1)
        if d == 1
            cp_dim = ones(Int64, ev[d]+1)
            push!(cp, cp_dim)
        else
            cp_dim = collect(StepRange(1, Int8(d), ev[d]*d+ 1))
            cp_dim = convert(Array{Int64,1}, cp_dim)
            push!(cp, cp_dim)
        end
    end
    return cp
end

function create_rv(
    ev,
    index2simplex,
    simplex2index
    )
    # create "rv" for CSC format input to Eirene
    
    """
    --- input ---
    ev: (array)
        ev[d]: number of (d-1) dimensional simplices
    --- output ---
    rv: (array) rv vector to be used as input to Eirene
    """
    maxdim = size(ev, 1)
    rv = []
    for d = 0:maxdim - 1
        # 0-simplices
        if d == 0
            push!(rv, Array{Int64,1}(undef, 0))

        # 1-simplices
        elseif d == 1
            rv_d = Array{Int64,1}(undef, 0)
            for i = 1:ev[d+1]
                append!(rv_d, index2simplex[(i,d)])
            end
            push!(rv, rv_d)

        # higher-dim simplices
        else
            rv_d = Array{Int64,1}(undef, 0)
            for i = 1:ev[d+1]
                boundary_idx = [simplex2index[item][1] for item in combinations(index2simplex[i,d], d)]
                append!(rv_d, boundary_idx)
            end
            push!(rv, rv_d)
        end
    end
    return rv
end

function create_CSC_input_from_fv(
    fv,
    index2simplex,
    simplex2index)
    
    # Find input for Eirene in CSC format (rv, cp, ev, fv), given the face values
    
    ##### ev: number of simplices in each dimension
    ev = [size(item,1) for item in fv]
    
    ##### boundary matrix CSC format: columns
    cp = create_CP(ev)
    
    ##### boundary matrix CSC format: rows
    rv = create_rv(ev, index2simplex, simplex2index)
    return rv, cp, ev, fv
    
end

function create_CSC_input_from_distance(
    n_simplices::Array{Int64,1}, 
    index2simplex::Dict, 
    simplex2index::Dict, 
    D_Y::Array{Float64,2})
    # Find input for Eirene in CSC format (rv, cp, ev, fv).
    # Use distance matrix D_Y to find the birth time of each simplex.
    """
    --- input ---
    n_simplices: output of function 'index_simplices_Zpsi'
    index2simplex: output of function 'index_simplices_Zpsi'
    simplex2index: output of function 'index_simplices_Zpsi'
    D_Y: distance matrix used for C_Y
    --- output ---
    rv: rv for Zpsi
            row values of the boundary matrix
    cp: cp for Zpsi
    ev: ev for Zpsi
            number of simplices in each dimension
    fv: fv for Zpsi
            birth time of simplices according to D_Y
    """

    ##### ev: number of simplices in each dimension
    maxdim = size(n_simplices,1)
    ev = n_simplices
    
    ##### boundary matrix CSC format: columns
    cp = create_CP(ev)
    
    ##### boundary matrix CSC format: rows
    rv = create_rv(ev, index2simplex, simplex2index)

    ##### fv: birth time of each simplex using D_X
    fv = []

    # 0-simplices have birth time of zero.
    fv_0 = zeros(ev[1])
    push!(fv, fv_0)

    # higher dimensional simplices
    for d = 1:maxdim-1
        fv_d = []
        for i = 1:ev[d+1]
            # find vertex representation of simplex using original indexing
            simplex = index2simplex[(i,d)]

            # find its birth time according to D_Y
            push!(fv_d, find_birth_time(simplex, D_Y))
        end

        push!(fv, fv_d)
    end
    return rv, cp, ev, fv
end

"""
    create_CSC_input_from_W_distance(;<keyword arguments>)
Given a complex, create Witness complex based on a given distance matrix. That is, use the distance matrix of the Witness complex to find the birth time. 
Return the input for Eirene in CSC format (rv, cp, ev, fv).

### Arguments
- `n_simplices`: number of simplices in complex
- `index2simplex` (dict): simplex index to simplex
- `simplex2index` (dict): simplex (as a list of vertices) to simplex
- `D_W`: the Witness distance matrix used to determine the birth time of each simplex

### Outputs
- `rv`: rv for W
    row values of the boundary matrix
- `cp`: cp for W
- `ev`: ev for W
    number of simplices in each dimension
- `fv`: fv for W
    birth time of simplices according to `D_W`
"""
function create_CSC_input_from_W_distance(
    n_simplices::Array{Int64,1}, 
    index2simplex::Dict, 
    simplex2index::Dict, 
    D_W::Array{Float64,2})

    ##### ev: number of simplices in each dimension
    maxdim = size(n_simplices,1)
    ev = n_simplices
    
    ##### boundary matrix CSC format: columns
    cp = ext.create_CP(ev)
    
    ##### boundary matrix CSC format: rows
    rv = ext.create_rv(ev, index2simplex, simplex2index)

    ##### fv: birth time of each simplex 
    fv = []

    # 0-simplices have birth time of zero.
    fv_0 = zeros(ev[1])
    push!(fv, fv_0)

    # higher dimensional simplices
    for d = 1:maxdim-1
        fv_d = []
        for i = 1:ev[d+1]
            # find vertex representation of simplex using original indexing
            simplex = index2simplex[(i,d)]

            # find its birth time according to D_W
            push!(fv_d, minimum(maximum(D_W[simplex,:], dims = 1)))
        end

        push!(fv, fv_d)
    end
    return rv, cp, ev, fv
end


#################################################################################
# Other helper functions
#################################################################################
function get_vertex_perm(C::Dict)
    # Eirene permutes the vertex indices.
    # Get vertex permutation information from Eirene.
    # note C["nvl2ovl"][i] = k, where i is the Eirene index and k is the original vertex index. 

    """
    --- input ---
    C: (dict) output from Eirene
    --- output ---
    v_perm: (arr) v_perm[k] = i, where k: original vertex index, 
            and i: corresponding Eirene index.
    
            ex) If C = Eirene(distance_matrix), 
            then k corresponds to the kth row/column of distance_matrix
    """

    n_vertices = size(C["nvl2ovl"],1)
    v_perm = zeros(Int64,n_vertices)
    for i=1:n_vertices
        idx = findall(x->x==i, C["nvl2ovl"])[1]
        v_perm[i] = idx
    end
    return v_perm
end


function select_simplices_at_psi(
    C::Dict, 
    psi::Float64;  
    maxdim::Int64 = 2)
    # Select the simplices in C that exist at parameter psi and return their indices. 
    # Note: Indexing is according to C, which is Eirene's internal indexing of simplices 

    """
    --- input ---
    C: output of running eirene
    psi: parameter
    maxdim: maximum dimension of simplices
    --- output ---
    simplices: array of simplices that exist at parameter delta 
            simplices[i] is an array of (i-1) dimensional simplices that exist at parameter delta
    """
    rv, cp, fv = Eirene_var.eirened2complex(C)
    simplices = []
    for d = 1:maxdim + 1
        push!(simplices, findall(x -> x <= psi, fv[d]))
    end
    return simplices
end

function create_idx_to_simplex_V(
    C::Dict; 
    maxdim::Int64 = 2)
    # Creates a dictionary from index of a simplex to its vertex representation, 
    # where the vertices are indexed according to C. 
    """
    --- input ---
    C: (Dict) output of Eirene. The maximum dimension of C must be at least "maxdim"-1 
    maxdim: (int) maximum dimension of simplices to consider
    --- output ---
    idx_smplx: (Dict) with (index, dimension) as key and its vertex representation as value
                ex: idx_smplx[(idx, dim)] = [v0, v1, ..., v_dim]

    --- note ---
    note: dim >= 1
    !!!!! NOTE !!!!! : The vertices are indexed according to Eirene's internal representation in C.
        (not the user-specified order)
    """
    # check that Eirene dictionary C has enough dimensions
    if C["input"]["maxdim"] < maxdim - 1
        throw(error("maxdim of Eirene output is too small."))
    end

    # get boundary matrix
    rv, cp, fv = Eirene_var.eirened2complex(C)
        
    idx_smplx = Dict()
    for d = 1:maxdim
        # 1-simplices
        if d == 1
            n = size(fv[d+1], 1)
            for i = 1:n
            # find vertex representation [v0, v1]
            idx_smplx[(i, 1)] = rv[2][i * 2 - 1:i * 2] 
            end
            
        # higher dimensional simplices
        else
            n = size(fv[d+1],1)
            for i = 1:n
                # find simplex as a list of its boundary
                # simplex = [b_1, b_2, ... ,b_{d+1}], where each b_i is the index 
                # of the i-th boundary simplex.
                simplex = rv[d+1][i * (d + 1) - d: i*(d+1)] 
                
                # find vertices of the simplex
                vertices = []
                # it suffices to just consider two of the boundary cells
                for item in simplex[1:2]
                    append!(vertices, idx_smplx[(item, d-1)])
                end
                idx_smplx[(i,d)] = unique(vertices)
            end
        end
    end
    return idx_smplx
end

function index_simplices_Zpsi(
    C_Z::Dict, 
    psi::Float64;
    maxdim::Int64 = 1)
    # Find simplices in the fixed complex Z_psi and devise an indexing scheme. 
    """
    --- input ---
    C_Z: Dictionary output of Eirene
    psi: parameter for building Z_psi
    maxdim: maximum dimension of homology of interest 

    --- output ---
    n_simplices: (array)
            n_simplices[i] is the number of (i-1) dimensional simplices of C_Z.
    Zpsi_index2simplex: (dict) given an index and dimension, returns the vertex representation of the simplex.
            Zpsi_index2simplex[(idx, dim)] = [v0, ... , vn]
    Zpsi_simplex2index: (dict) reverse of Zpsi_index2simplex
            Zpsi_simplex2index[[v0, ... , vn]] = idx

    !!!!! NOTE !!!!!
    In Zpsi_index2simplex and Zpsi_simplex2index, the vertices are ordered according to the original indexing, and not C_Z.
    That is, if C_Z was obtained from distance matrix D_Z, vertex "i" corresponds to the i-th row and column of D_Z
    """

    ##### 1. Find all simplices of C_Z, their index (according to C_Z), and their vertex representation (using original indexing)
    Z_index2simplexV_eirene = create_idx_to_simplex_V(C_Z; maxdim = maxdim + 1)

    # express each simplex in C_Z using original vertices
    v_perm = C_Z["nvl2ovl"]
    Z_index2simplexV_orig = Dict(key => v_perm[value] for (key, value) in Z_index2simplexV_eirene)

    # reverse
    Z_simplex2indexV_orig = Dict(val => key for (key, val) in Z_index2simplexV_orig)

    ##### 2. create list of simplices that exist in Zpsi
    # note: "simplices_psi" is a list of indices, where the index corresponds to Eirene's indexing of C_Z
    simplices_psi = select_simplices_at_psi(C_Z, psi, maxdim = maxdim + 1)

    ##### 3. Create indexing of simplices in Z_psi

    # 0-simplices
    n0 = size(simplices_psi[1],1)
    Zpsi_simplex2index = Dict( [i] => v_perm[i] for i =1:n0)

    # higher dimensional simplices
    for d = 1:size(simplices_psi, 1)-1
        for (index, value) in enumerate(simplices_psi[d+1])
            simplex = sort(Z_index2simplexV_orig[(value, d)])
            Zpsi_simplex2index[simplex] = index
        end
    end

    # reverse dictionary
    Zpsi_index2simplex = Dict((val, size(key,1)-1) => key for (key, val) in Zpsi_simplex2index)

    # find number of simplices in each dimension
    n_simplices = [size(item,1) for item in simplices_psi]

    return n_simplices, Zpsi_index2simplex, Zpsi_simplex2index
end


function find_birth_time(
    simplex::Array, 
    D::Array{Float64,2})
    # Given a simplex (written as a list of vertices), find its birth time using distnace matrix D
    """
    --- input ---
    simplex: array of vertices.
            ex) simplex = [1,2,3]
    D: distance matrix
    --- output ---
    time: (float) birth time of given simplex. 
    """
    time = 0
    for item in combinations(simplex, 2)
        time = max(time, D[item[1], item[2]])
    end
    return time
end





function find_linear_combinations(
    candidates)
    # Given an array of integers (corresponding to vectors), find all non-trivial linear combinations
    """
    --- input ---
    candidates: (array) of integers
    --- output ---
    linear_comb: (array) of possible non-trivial linear combinations of candidates
    """

    linear_comb = []
    for i = 1:size(candidates, 1)
        comb = collect(combinations(candidates, i))
        append!(linear_comb, comb)
    end
    return linear_comb
end

function select_odd_count(
    orig_list::Array)
    # given an array, return the elements that occur odd number of times.
    """
    --- input ---
    orig_list: (N-element array)
    --- output ---
    new_list: (M-element array) containing only the elements that occur odd number of times
            in orig_list.
    """

    count = countmap(orig_list)
    new_list = [item for item in orig_list if count[item] % 2 != 0]
    return unique(new_list)
end

function express_x_CZpsi(
    C_Z::Dict, 
    class_num::Int64, 
    Z_psi_simplex2index::Dict;
    dim::Int64=1)
    # Express [x] in H_n(C_Z) as a chain (list of indices) using the indexing of C_Z_psi)
    """
    --- input ---
    C_Z: (dict) output of Eirene
    class_num: (int) class number of cycle [x] of interest
    Z_psi_simplex2index: (dict) simplex indexing of Z_psi.
                            output of function 'index_simplices_Zpsi'
    dim: (int) dimension of cycle [x]
    --- output ---
    """
    # express [x] using vertices.
    # note: cycle_rep format is "vertex x simplex", where vertices are indexed using the original indexing
    cycle_rep = classrep(C_Z, class = class_num, dim = dim)

    chain = []
    #  for each simplex in [x]
    for i = 1:size(cycle_rep, 2)
        # express simplex using original vertices
        simplex = cycle_rep[:,i]

        # find index in C_Z_filt
        idx = Z_psi_simplex2index[sort(simplex)]
        append!(chain, idx)
    end    
    return chain
end
    

function bounding_chain(C;chain=zeros(Int64,0),dim=1)
    # Check if a given chain is a boundary
    
    ##### CREDITS: This function was written by Greg Henselman-Petrusek. To be included in a future version of Eirene #####
    ##### https://github.com/Eetion/Eirene.jl #####
    
	if isempty(chain)
		return zeros(Int64,0)
	elseif !isempty(Eirene_var.chainboundary(C,chain=chain,dim=dim))
		print("there is no bounding chain"); return nothing
	else
		sd 			= 	dim+1;
		Lrv 		= 	C["Lrv"][sd+1]
		Lcp 		= 	C["Lcp"][sd+1]
		Rrv 		= 	C["Rrv"][sd+1]
		Rcp 		= 	C["Rcp"][sd+1]
		numrows 	= 	Eirene_var.boundarycorank(C,dim=dim)
		translator 	= 	zeros(Int64,Eirene_var.complexrank(C,dim=dim))
		translator[C["tid"][sd+1]] = 1:numrows
		rv 			= 	findall(!iszero,translator[chain])
		rv 			=  	translator[chain[rv]]
		cp 			= 	[1,length(rv)+1]
		rv,cp 		=  	Eirene_var.spmmF2silentLeft(Lrv,Lcp,rv,cp,numrows)
		rv,cp 		=  	Eirene_var.spmmF2silentLeft(Rrv,Rcp,rv,cp,numrows)
		#
		# recall that plo = the (ordered) subvector consisting of
		# the first (rank(boundary operator)) elements of tid
		#
		if maximum(rv) > length(C["plo"][sd+1])
			#print("there is no bounding chain");  # NOTE: Iris removed print
			return nothing
		else
			ocg2rad = 	C["ocg2rad"]
			grain 	= 	C["grain"][sd+1]
			names 	= 	C["phi"][sd+1]
			return 	names[rv]
		end
	end
end





function simplex_to_index(simplex, C)
    # Given a n-simplex [u_0, u_1, ... , u_n], where each u_i is Eirene-indexed,
    # find the index of simplex in C_n according to eirene output C.
    
    """
    --- input ---
    simplex (arr) : [u_0, u_1, ... , u_k], where each u_i is indexed 
                according to Eirene.
                u_0 < u_1 < ... < u_k.
    C (dict): output from Eirene
    --- output ---
    index (int): index of given simplex in C_n according to C.

    --- example use ---
    # 1. get permutation info
    v_perm = get_vertex_perm(C)
    # 2. rewrite simplex in terms of 0-simplices (indexed according to eirene)
    simplex = [96, 183, 188]
    simplex_eirene = sort([v_perm[i] for i in simplex])
    # 3. find simplex index
    simplex_to_index(simplex_eirene, C)
    """ 
    dim = size(simplex,1) - 1

    # given simplex_eirene =[u_0, u_1, ... , u_n], find the portion of row-indices for column u0
    rowvals = Eirene_var.crows(C["firstv"][dim+1], C["farfaces"][dim+1], simplex[1])

    # find location of [u_1, ..., u_n] in rowvals.
    # that is, find the index of [u_1, ..., u_n] among C_{n-1}.
    
    # if simplex is 1-dimensional
    if dim == 1
        rowloc = findall(x->x==simplex[2], rowvals)[1]
        
    # if simplex is higher dimensional
    else
        rowloc0 = simplex_to_index(simplex[2:end], C)
        rowloc = findall(x->x==rowloc0, rowvals)[1]
    end

    # index of simplex according to Eirene
    index = C["firstv"][dim+1][simplex[1]] + rowloc - 1
    
    return index
end



function chain_to_index(
    chain, 
    C::Dict)
    # Given an n-dimensional chain expressed as a list of simplices 
    # chain = [simplex_1, simplex_2, ... , simplex_k],
    # where each simplex_j = [v_0, ... ,v_n], a list of its (n+1) vertices,
    # rewrite the chain as a list of integers [i_1, ..., i_k], 
    # where each $i_j$ is the index of simplex_j in C_n according to Eirene. 

    """
    --- input ---
    chain (arr) : [s1, s2, ... , sn], 
                where each s_j is a n-simplex of form [v_0, v_1, ..., v_n],
                a list of its (n+1) vertices
    C (dict): output from Eirene

    --- output ---
    chain_idx (arr): [i_1, ... ,i_n], where each i_j is the index of s_j according to C

    --- example ---
    test_chain =[[96, 183, 188],[99, 111, 188]]
    chain_idx = chain_to_index(test_chain, C)
    print(chain_idx)

    # check
    dim = size(test_chain[1],1)-1
    for item in chain_idx
        simplex = EireneVar.incidentverts(C["farfaces"], C["firstv"], dim+1, [item])
        simplex = sort(C["nvl2ovl"][simplex])
        print(simplex)
    end
    """
    # get permutation of vertices
    v_perm = get_vertex_perm(C)

    chain_idx = []
    for simplex in chain
        simplex_eirene = sort([v_perm[i] for i in simplex])
        try 
            simplex_idx = simplex_to_index(simplex_eirene, C)
            append!(chain_idx, simplex_idx)
        catch
            print("chain doesn't exist in eirene output")
            print(simplex)
            return nothing
        end
    end

    return chain_idx
end


"""
    compute_distance(P, Q)
Compute the Euclidean distance matrices

### Arguments
- `P`: (array) of points in one region. Rows: points, columns: coordinates.
- `Q`: (array) of points in second region. Rows: points, columns: coordinates.

### Outputs
- `D_P`: (array) distance matrix among points in P
- `D_Q`: (array) distance matrix among points in Q
- `D_P_Q`: (array) distance matrix among points in P (row) and Q (column)
- `D_Q_P`: (array) distance matrix among points in Q (row) and P (column)
"""
function compute_distance(P, Q)
    
    # number of points in P
    n_P = size(P,1)
    
    # gather all points 
    X = vcat(P, Q)

    # compute distance
    D = Distances.pairwise(Euclidean(), X, X, dims = 1)

    # Define submatrices 
    D_P = D[1:n_P, 1:n_P]
    D_Q = D[n_P+1:end, n_P+1:end]
    D_P_Q = D[1:n_P, n_P+1:end]
    D_Q_P = D[n_P+1:end, 1:n_P]
    
    return D_P, D_Q, D_P_Q, D_Q_P
end


function plot_PD(barcode; 
    highlight = [], highlight_color = :deeppink2, cutoff = nothing, 
    pd_min = nothing, pd_max = nothing, threshold_lw = 5, diagonal_lw = 5,
    kwargs...)
points = barcode


if size(barcode,1) == 0
    # plot diagonal line
    p = plot([0, 1], [0, 1], 
    labels ="", 
    linewidth = diagonal_lw,
    framestyle = :box,
    xlims = (0,1),
    ylims = (0,1),
    aspect_ratio = :equal,
    color = "grey"; 
    kwargs...)
    return p
end
# find index of points with death parameter == death
idx = findall(x -> x == Inf, points[:,2])

# plot points with death < Inf
idx2 = [i for i in 1:size(points,1) if i ∉ idx]
p = scatter(points[idx2,1], points[idx2,2]; kwargs..., color = "grey", labels = "", alpha = 0.5, hover = idx2)

# find max death value
max_death = maximum(points[idx2, 2])

# plot points with death parameter == Inf
death_plot = ones(size(idx,1)) * max_death

scatter!(p, points[idx,1], death_plot, aspect_ratio = :equal, legend=:bottomright, labels="", color ="red"; kwargs...)

# plot diagonal line
if pd_max == nothing
    
    min_birth = minimum(points[:,1]) * 0.8
    max_death = max_death * 1.1
    plot!(p, [min_birth, max_death], [min_birth, max_death], 
        labels ="", 
        linewidth = diagonal_lw,
        framestyle = :box,
        xlims = (min_birth, max_death),
        ylims = (min_birth, max_death),
        color = "grey"; 
        kwargs...)
else
    max_death = pd_max
    min_birth = pd_min
    plot!(p, [min_birth, max_death], [min_birth, max_death], 
        labels ="", 
        linewidth = diagonal_lw,
        framestyle = :box,
        xlims = (min_birth, max_death),
        ylims = (min_birth, max_death),
        color = "grey"; 
        kwargs...)
end

 # if highlight is provided, color specific points with the given color
if highlight != []
     scatter!(points[highlight,1], points[highlight,2]; kwargs..., color = highlight_color, labels = "", hover = highlight)
end

# plot the cutoff (dotted line) if provided
if cutoff != nothing
    f(x) = x + cutoff
    plot!(f, linestyle = :dash, c = "black", label = "", hover = false, linewidth = threshold_lw)
end

return p
end


function get_1simplices(D_param)
    n_cols = size(D_param, 2)
    one_simplices = []
    for i = 1:n_cols
        rows = findall(x -> x == 1, D_param[:,i])
        append!(one_simplices, combinations(rows, 2))
    end
    
    return one_simplices
end

function get_2simplices(D_param)
    n_cols = size(D_param,2)
    two_simplices = []
    for i = 1:n_cols
        ones = findall(D_param[:,i])
        append!(two_simplices, collect(combinations(ones,3)))
    end
    
    return unique(two_simplices)
end

"""
    plot_Dowker_complex
Given a dissimilarity matrix D, compute the Dowker complex at parameter `param`.
Uses the rows as potential vertex set. Must provide `PC`, the locations of vertex set corresponding to the rows.

### Inputs
- `D`: Dissimilarity matrix
- `param`: Parameter. Builds the Dowker complex at this parameter.
- `PC`: An array of size (n,2), where n is the number of rows of D. 
        Provides (x,y)-coordinates of vertices corresponding to the  rows of D. 

### Outputs
- a plot
"""
function plot_Dowker_complex(D, param, PC; 
                             show_2simplex = false, 
                             show_unborn_vertices = false,
                             c = "#29A0B1",
                            kwargs...)
    n_rows, n_cols = size(D)
    D_param = D .< param
    
    p = plot()
    
    # plot 2-simplices
    if show_2simplex == true
        two_simplices = get_2simplices(D_param)
        for simplex in two_simplices
            plot!(p, Plots.Shape([(PC[i,1], PC[i,2]) for i in simplex]), 
                 label="", legend = :false, c = c, alpha = 0.1
                 )
        end
    end
    
    # plot 1-simplices
    one_simplices = get_1simplices(D_param)
    for simplex in one_simplices
        plot!(p, [PC[simplex[1],1], PC[simplex[2],1]], [PC[simplex[1],2], PC[simplex[2],2]], label = "", c = :black) 
    end
    
    # plot 0-simplices
    idx0 = findall(x -> x != 0, vec(sum(D_param, dims = 2)))
    scatter!(p, PC[idx0,1], PC[idx0, 2]; label = "", frame = :box, ticks = [], c = c, aspect_ratio = :equal, kwargs...)
    
    # plot unborn vertices 
    if show_unborn_vertices == true
        idx_unborn = findall(x -> x == 0, vec(sum(D_param, dims = 2)))
        scatter!(p, PC[idx_unborn,1], PC[idx_unborn, 2], label = "",
                 markershape = :xcross,
                 markerstrokewidth = 3,
                 c = c) 
    end
    return p
end

end
