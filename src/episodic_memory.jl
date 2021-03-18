struct AppearanceFeature end
struct BothFeatures end
mutable struct SnapshotContent{T}
    weight::Float64 # interpretation: weight from feature to food type
    t0::T           # used to compute the corresponding age bucket.
end
abstract type AbstractSnapshots end
struct SnapshotSimple{F, K, T} <: AbstractSnapshots
    mapping::Dict{K, Dict{Food, SnapshotContent{T}}}
end
function SnapshotSimple(T; feature = AppearanceFeature, kwargs...)
    K = feature == BothFeatures ? Tuple{Int, Int} : Int
    SnapshotSimple{feature, K, T}(Dict{K, Dict{Food, SnapshotContent{T}}}())
end
_inspect_features(s::SnapshotSimple{AppearanceFeature}, tray) = tray.appearance
_inspect_features(s::SnapshotSimple{BothFeatures}, tray) = (tray.position, tray.appearance)
cache_features(::Any, tray) = tray.position
function isdepleted(m, f, t = nothing)
    if t === nothing
        for content in values(m[f])
            content.weight > 0 && return false
        end
    else
        m[f][t].weight > 0 && return false
    end
    true
end
function Base.in(s::SnapshotSimple, tray::Tray, typ = nothing)
    f = _inspect_features(s, tray)
    haskey(s.mapping, f) && (typ === nothing || haskey(s.mapping[f], typ)) && !isdepleted(s.mapping, f, typ)
end
function foodtypes(s::SnapshotSimple, tray)
    f = _inspect_features(s, tray)
    foodtypes = Food[]
    !haskey(s.mapping, f) && return foodtypes
    for (k, v) in s.mapping[f]
        v.weight <= 0 && continue
        push!(foodtypes, k)
    end
    foodtypes
end
function add_snapshot_item!(s::SnapshotSimple, tray, typ, t)
    f = _inspect_features(s, tray)
    if !in(s, tray)
        s.mapping[f] = Dict()
    end
    if haskey(s.mapping[f], typ)
        s.mapping[f][typ].weight += 1 # rand() # fuzzy consolidation
        s.mapping[f][typ].t0 = t
    else
        s.mapping[f][typ] = SnapshotContent(1., t)
    end
    s
end
function remove_snapshot_item!(s::SnapshotSimple, tray, typ)
    if in(s, tray)
        f = _inspect_features(s, tray)
		if haskey(s.mapping[f], typ)
            s.mapping[f][typ].weight -= 1 # rand() # fuzzy deconsolidation
		end
	end
    s
end
adjust_snapshot_features!(s, tray) = nothing
memoryage(s::SnapshotSimple, tray, typ) = s.mapping[_inspect_features(s, tray)][typ].t0
function memory_count(s::SnapshotSimple, tray, typ)
    n = 0.
    if in(s, tray, typ)
        n = s.mapping[_inspect_features(s, tray)][typ].weight
    end
    n
end

struct Snapshots{T} <: AbstractSnapshots
    weights::Vector{Float64}
    snapshots::T
    alpha_adjust::Float64
end
function Snapshots(T; kwargs...)
    p = paramvals(Snapshots; kwargs...)
    Snapshots([p.positionweight, p.appearanceweight],
              SnapshotSimple(T, feature = :both),
              p.alpha_adjust)
end
function Base.in(s::Snapshots, tray, typ = nothing)
    f = _inspect_features(s.snapshots, tray)
    if haskey(s.snapshots.mapping, f) &&
        (typ === nothing || haskey(s.snapshots.mapping[f], typ))
        return !isdepleted(s.snapshots.mapping, f, typ)
    end
    # fallback
    for k in keys(s.snapshots.mapping)
        tray.position == k[1] || tray.appearance == k[2] && return true
    end
    false
end
function params(::Type{Snapshots})
    (positionweight = (searchrange = (.1, .9), default = .6),
     appearanceweight = (searchrange = (.1, .9), default = .4),
     alpha_adjust = (searchrange = (.8, 1.), default = .9))
end
function _recalled_tray(s, tray)
    length(s.snapshots.mapping) == 0 && return nothing
    positions = Tuple{Int, Int}[]
    appearances = Tuple{Int, Int}[]
    for k in keys(s.snapshots.mapping)
        if tray.position == k[1]
            push!(positions, k)
        elseif tray.appearance == k[2]
            push!(appearances, k)
        end
    end
    if length(positions) > 0
        if length(appearances) > 0
            return Tray(rand(wsample([positions, appearances], s.weights))...)
        else
            return Tray(rand(positions)...)
        end
    elseif length(appearances) > 0
        return Tray(rand(appearances)...)
    end
    nothing
end
function foodtypes(s::Snapshots, tray)
    if haskey(s.snapshots.mapping, _inspect_features(s.snapshots, tray))
        foodtypes(s.snapshots, tray)
    else
        r = _recalled_tray(s, tray)
        r === nothing && return Food[]
        foodtypes(s.snapshots, r)
    end
end
function adjust_snapshot_features!(s::Snapshots, tray)
    for k in keys(s.snapshots.mapping)
        if tray.position == k[1]
            s.weights[1] *= s.alpha_adjust
        elseif tray.appearance == k[2]
            s.weights[2] *= s.alpha_adjust
        end
    end
end
function remove_snapshot_item!(s::Snapshots, tray, typ)
    if haskey(s.snapshots.mapping, _inspect_features(s.snapshots, tray))
        remove_snapshot_item!(s.snapshots, tray, typ)
    else
        r = _recalled_tray(s, tray)
        r === nothing && return nothing
        remove_snapshot_item!(s.snapshots, r, typ)
    end
end
function add_snapshot_item!(s::Snapshots, tray, typ, t)
    add_snapshot_item!(s.snapshots, tray, typ, t)
end
function cache_features(s::Snapshots, tray)
    cache_features(s.snapshots, tray)
end
function memoryage(s::Snapshots, tray, typ)
    if haskey(s.snapshots.mapping, _inspect_features(s.snapshots, tray))
        memoryage(s.snapshots, tray, typ)
    else
        r = _recalled_tray(s, tray)
        memoryage(s.snapshots, r, typ)
    end
end
function memory_count(s::Snapshots, tray, typ)
    if haskey(s.snapshots.mapping, _inspect_features(s.snapshots, tray))
        memory_count(s.snapshots, tray, typ)
    else
        r = _recalled_tray(s, tray)
        r === nothing && return 0
        memory_count(s.snapshots, r, typ)
    end
end

# api
# function Base.in(snapshots::AbstractSnapshots, tray::Tray, typ = nothing) end
# function foodtypes(snapshots, tray) end
# function remove_snapshot_item!(snapshots, tray, typ) end
# function add_snapshot_item!(snapshots, tray, typ, t) end
# function _inspect_features(snapshots, tray) end
# function cache_features(snapshots, tray) end
# function memoryage(snapshots, tray, typ) end
# function memory_count(snapshots, tray, typ) end
# adjust_snapshot_features!(snapshots, tray) end
