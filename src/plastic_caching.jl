#TODO: unify datastructure for trayweights and snapshots
struct MultiplicativeFreshness end
struct AdditiveFreshness end
struct StaticCachingAgentParams{F}
	inspectbias::Float64
    inspectslope::Float64
end
freshnessweight(::StaticCachingAgentParams{AdditiveFreshness}, h, f) = h + f - 1
freshnessweight(::StaticCachingAgentParams{MultiplicativeFreshness}, h, f) = h * f
inspectbias(s::StaticCachingAgentParams) = s.inspectbias
inspectslope(s::StaticCachingAgentParams) = s.inspectslope
hastrayweights(::Any) = false
struct AbundanceDecrease
    rate::Float64
end
AbundanceDecrease(p::NamedTuple) = AbundanceDecrease(p.abundancedecrease)
struct HungerIncrease{T}
    hungerincreasethreshold::Float64
    hungerincreasetimeconstant::T
end
struct PlasticCachingAgentParams{S,T,C}
    staticparams::S
    maxcacheweight::Float64
	pilferdecrease::Float64
	retrievalrewardfactor::Float64
    degradedecrease::Float64
    dontcare::Float64
	compensatorycaching::C
    fillleveldecreasemodel::T
	initialweight::Float64
    trayweights::Dict{Tuple{Int, Food},Float64}
end
freshnessweight(p::PlasticCachingAgentParams, h, f) = freshnessweight(p.staticparams, h, f)
inspectbias(s::PlasticCachingAgentParams) = inspectbias(s.staticparams)
inspectslope(s::PlasticCachingAgentParams) = inspectslope(s.staticparams)
hastrayweights(::PlasticCachingAgentParams) = true

struct SigmaFreshnessEstimator{T}
    freshnessparameter::Array{Float64, 1}
    memorytimeconstant::T
    freshness_sigma::Float64
    freshness_alpha::Float64
end
function SigmaFreshnessEstimator(; freshnessparameter = 10*ones(N_FOODTYPE),
                                 memorytimeconstant = 1u"d",
                                 freshness_alpha = 1.5,
                                 freshness_sigma = .1, kwargs...)
    SigmaFreshnessEstimator(freshnessparameter, memorytimeconstant,
                            freshness_sigma, freshness_alpha)
end
struct PiecewiseConstant{T}
    boundaries::Array{T, 1}
    values::Array{Float64, 1}
end
struct PiecewiseConstantFreshnessEstimator{T}
    freshness::Array{PiecewiseConstant{T},1}
    freshness_alpha::Float64
end
function show(io::IO, mime::MIME"text/plain", pcfe::PiecewiseConstantFreshnessEstimator)
    println(io, "PiecewiseConstantFreshnessEstimator")
    println(io, "learningrate: $(pcfe.freshness_alpha)")
    for f in instances(Food)
        println(io, "$f")
        println(io, "boundaries $((x -> x.val).(pcfe.freshness[Int(f)].boundaries))")
        println(io, "values     $(pcfe.freshness[Int(f)].values)")
    end
end
mutable struct PlasticCachingAgent{P,TF,Ts,T,Tsnap,Th} <: SpecSatAgents
    t::T
    params::Params
    specsatparams::Ts
    cacheparams::P
    freshnessestimator::TF
    snapshots::Tsnap
    hungermodel::Th
end
hastrayweights(agent::PlasticCachingAgent) = hastrayweights(agent.cacheparams)

function expectedfreshness(estimator::SigmaFreshnessEstimator, Δt, foodtypeindex)
    θ = estimator.freshnessparameter[foodtypeindex]
    δ = θ - Δt/estimator.memorytimeconstant
    (delta = δ, freshness = normcdf(δ/(θ * estimator.freshness_sigma)))
end
function expectedfreshness(agent, tray, typ)
    expectedfreshness(agent.freshnessestimator,
                      agent.t - memoryage(agent.snapshots, tray, typ),
                      foodindex(typ))
end
function update_freshnessparameter!(agent::PlasticCachingAgent{<:Any, <: SigmaFreshnessEstimator},
                                    tray, typ, freshness)
    δ, expected = expectedfreshness(agent, tray, typ)
    agent.freshnessestimator.freshnessparameter[foodindex(typ)] -=
        agent.freshnessestimator.freshness_alpha * abs(freshness - expected) * δ
end
function expectedfreshness(estimator::PiecewiseConstantFreshnessEstimator,
                           Δt, foodtypeindex)
    i = 1
    while estimator.freshness[foodtypeindex].boundaries[i] < Δt
        i += 1
    end
    estimator.freshness[foodtypeindex].values[i], i
end
function update_freshnessparameter!(agent::PlasticCachingAgent{<:Any, <:PiecewiseConstantFreshnessEstimator},
                                    tray, typ, freshness)
    expected, index = expectedfreshness(agent, tray, typ)
    f = agent.freshnessestimator
    f.freshness[foodindex(typ)].values[index] += f.freshness_alpha * (freshness - expected)
end
inspectweight(agent::PlasticCachingAgent, tray::Tray) = _inspectweight(agent, tray)
function _inspectweight(agent, tray)
    if !tray.closed
        weight = 0.
        s = inspectslope(agent.cacheparams)
        if in(agent.snapshots, tray)
            for foodtype in foodtypes(agent.snapshots, tray)
                if agent.t - memoryage(agent.snapshots, tray, foodtype) > 1u"hr"
                    j = foodindex(foodtype)
                    f, i = expectedfreshness(agent, tray, foodtype)
                    weight = max(weight,
                                 freshnessweight(agent.cacheparams,
                                                 s * geteatweight(agent, j),
                                                 f))
                else
                    return 0.
                end
            end
        end
        return weight + inspectbias(agent.cacheparams)
    end
    return 0.
end

function inspectcallback!(agent::PlasticCachingAgent, tray::Tray)
    # pilfer
    if in(agent.snapshots, tray)
        for foodtype in foodtypes(agent.snapshots, tray)
            id = (cache_features(agent.snapshots, tray), foodtype)
            if hastrayweights(agent) && haskey(agent.cacheparams.trayweights, id)
                agent.cacheparams.trayweights[id] *= agent.cacheparams.pilferdecrease
            end
            hasmemory(agent) && update_model!(agent, id[1], :pilfered)
            remove_snapshot_item!(agent.snapshots, tray, foodtype)
        end
        adjust_snapshot_features!(agent.snapshots, tray)
    end
end

function grow_weight(agent, cf, typ, γ)
    agent.cacheparams.trayweights[(cf, typ)] *= γ
    agent.cacheparams.trayweights[(cf, typ)] += (1 - γ) * agent.cacheparams.maxcacheweight
end
function inspectcallback!(agent::PlasticCachingAgent, tray::Tray, item)
	# degrade & hunger dependence
    typ = item.id
    if in(agent.snapshots, tray, typ)
        update_freshnessparameter!(agent, tray, typ, item.freshness)
        cf = cache_features(agent.snapshots, tray)
        if hastrayweights(agent) && haskey(agent.cacheparams.trayweights, (cf, typ))
            i = foodindex(typ)
            eatweight = geteatweight(agent, i)
            if item.freshness > .5 && eatweight > agent.cacheparams.dontcare
                γ = 1 - (1 - agent.cacheparams.retrievalrewardfactor) * eatweight
                grow_weight(agent, cf, typ, γ)
            else
                γ = agent.cacheparams.degradedecrease
                agent.cacheparams.trayweights[(cf, typ)] *= γ
            end
        end
        if hasmemory(agent)
            if item.freshness > .5
                update_model!(agent, cf, :rewarded)
            else
                update_model!(agent, cf, :degraded)
            end
        end
        remove_snapshot_item!(agent.snapshots, tray, typ)
    end
end

inittrayweights!(::Any, ::Any, ::Any) = nothing
function inittrayweights!(agent::PlasticCachingAgent{<:PlasticCachingAgentParams}, position, typ)
    if hastrayweights(agent) && !haskey(agent.cacheparams.trayweights, (position, typ))
            agent.cacheparams.trayweights[(position, typ)] =
                agent.cacheparams.initialweight
    end
    nothing
end

function cachecallback!(agent::PlasticCachingAgent, item, tray)
    typ = item.id
	add!(tray, typ, 1)
    modifycachemotivation!(agent.hungermodel, agent.specsatparams, typ)
    add_snapshot_item!(agent.snapshots, tray, typ, agent.t)
    inittrayweights!(agent, cache_features(agent.snapshots, tray), typ)
end

function abundancedecrease!(agent::PlasticCachingAgent{<:PlasticCachingAgentParams{<:Any, <:AbundanceDecrease}}, typ, trays)
	for tray in trays
        if !in(agent.snapshots, tray, typ)
            inittrayweights!(agent, cache_features(agent.snapshots, tray), typ)
            agent.cacheparams.trayweights[(cache_features(agent.snapshots, tray), typ)] *=
                    agent.cacheparams.compensatorycaching.rate
        end
	end
    nothing
end

struct FillLevelDecrease
    f::Float64
end
struct NoFillLevelDecrease end
cacheweight_decrease(::NoFillLevelDecrease, agent, tray, typ) = 0.
function cacheweight_decrease(f::FillLevelDecrease, agent, tray, typ)
    n = memory_count(agent.snapshots, tray, typ)
    f.f * (1 - f.f^n)
end
function cachefromtoweight(agent::PlasticCachingAgent, object::FoodItem, tray::Tray)
    typ = object.id
    inittrayweights!(agent, cache_features(agent.snapshots, tray), typ)
    j = foodindex(typ)
    witem = getcacheweight(agent, j)
    wtray = 0.
    if hastrayweights(agent)
        wtray = agent.cacheparams.trayweights[(cache_features(agent.snapshots, tray), typ)]
    end
    wtray + witem
end

function updatecacheweights!(agent::PlasticCachingAgent{<:PlasticCachingAgentParams{<:Any, <:Any, <:HungerIncrease}}, cage, Δt)
    length(cage.trays) == 0 && return nothing
    hungermodel = agent.hungermodel
    idxs = eachindex(hungermodel)
    hθ = agent.cacheparams.compensatorycaching.hungerincreasethreshold
    hτ = agent.cacheparams.compensatorycaching.hungerincreasetimeconstant
    ds = [duration_above_threshold(hungermodel, Δt, hθ, cage.ismdpresent, i)
          for i in idxs]
    for i in idxs
        d = ds[i]
        typ = agent.specsatparams.foodtypes[i]
        for tray in cage.trays
            iszero(d) && continue
            cf = cache_features(agent.snapshots, tray)
            inittrayweights!(agent, cf, typ)
            γ = exp(-d/hτ)
            grow_weight(agent, cf, typ, γ)
        end
    end
end
updatecacheweights!(::Any, ::Any, ::Any) = nothing
function updateagentstate!(agent::PlasticCachingAgent, cage, Δt)
    agent.t += Δt
    updatecacheweights!(agent, cage, Δt)
    updatehunger!(agent, cage, Δt)
end
