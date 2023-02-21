function foodlookup(foodtypes)
    lookup = zeros(Int, N_FOODTYPE)
    for (i, v) in pairs(foodtypes)
        lookup[foodindex(v)] = i
    end
    tuple(lookup...)
end

struct HungerModulatedCaching end
update_cachemod!(::HungerModulatedCaching, ::Any, ::Any, ::Any) = nothing
modifycachemotivation!(::HungerModulatedCaching, ::Any, ::Any) = nothing
struct CacheModulatedCaching
    cachedecreasescale::Float64
    cachemotivation::Vector{Float64}
end
function update_cachemod!(c::CacheModulatedCaching, Δt, τ, i)
    c.cachemotivation[i] = 1 + (c.cachemotivation[i] - 1) * exp(-Δt/τ)
end
function modifycachemotivation!(c::CacheModulatedCaching, specsat, typ)
    i = lookup(specsat, typ)
    c.cachemotivation[i] -= c.cachemotivation[i] * c.cachedecreasescale
end
struct CacheModulatedCaching2{T,N}
    hungertimeconstant::T
    digestiontimeconstant::T
    digestionduration::T
    update_value::NTuple{N,Float64}
    hunger::Vector{Float64}
    stomach::Vector{Float64}
end
function update_cachemod!(c::CacheModulatedCaching2, Δt, ::Any, i)
    update!(c, Δt, false, i)
end
function modifycachemotivation!(c::CacheModulatedCaching2, specsat, typ)
    i = lookup(specsat, typ)
    c.stomach[i] += c.update_value[i]
end

struct Hunger{C,T}
    hungertimeconstant::T
    digestiontimeconstant::T
    digestionduration::T
    hunger::Vector{Float64}
    stomach::Vector{Float64}
    cachemodulation::C
end

hunger(h::Hunger) = h.hunger
Base.eachindex(h::Hunger) = eachindex(h.hunger)
hunger(agent) = hunger(agent.hungermodel)
eat!(h::Hunger, nutritionvalue::Number, i) = h.stomach[i] += nutritionvalue
modifycachemotivation!(h::Hunger, specsat, typ) = modifycachemotivation!(h.cachemodulation, specsat, typ)
function eat!(h, nutritionvalue)
    for i in eachindex(nutritionvalue)
        eat!(h, nutritionvalue[i], i)
    end
end
function duration_above_threshold(h::Hunger, Δt, θ, ismdpresent, i)
    if ismdpresent
        return duration_above_threshold_decr(h.hunger[i], θ, Δt, h.digestiontimeconstant)
    end
    if h.stomach[i] > 0.
        timeuntildigested = h.stomach[i] * h.digestionduration
        if Δt > timeuntildigested
            d_above = duration_above_threshold_decr(h.hunger[i], θ,
                                                    timeuntildigested,
                                                    h.digestiontimeconstant)
            hungerdec = h.hunger[i] * exp(-timeuntildigested/h.digestiontimeconstant)
            newΔt = Δt - timeuntildigested
            d_above + duration_above_threshold_incr(hungerdec, θ, newΔt, h.hungertimeconstant)
        else
            duration_above_threshold_decr(h.hunger[i], θ, Δt, h.digestiontimeconstant)
        end
    else
        duration_above_threshold_incr(h.hunger[i], θ, Δt, h.hungertimeconstant)
    end
end
function duration_above_threshold_decr(h, θ, Δt, τ)
    h ≤ θ && return 0.0u"minute"
    return min(τ * log(h/θ), Δt)
end
function duration_above_threshold_incr(h, θ, Δt, τ)
    h ≥ θ && return Δt
    return max(0.0u"minute", Δt + τ * log((θ-1)/(h-1)))
end
function update!(h, Δt, ismdpresent = false)
    for i in eachindex(h)
        update!(h, Δt, ismdpresent, i)
        update_cachemod!(h.cachemodulation, Δt, h.hungertimeconstant, i)
    end
    hunger(h)
end
decreasehunger!(h, Δt, i) = h.hunger[i] *= exp(-Δt/h.digestiontimeconstant)
function increasehunger!(h, Δt, ismdpresent, i)
    ismdpresent && return decreasehunger!(h, Δt, i)
    h.hunger[i] = 1 + (h.hunger[i] - 1) * exp(-Δt/h.hungertimeconstant)
end
function update!(h::Union{Hunger, CacheModulatedCaching2}, Δt, ismdpresent, i)
    if h.stomach[i] > 0.
        timeuntildigested = h.stomach[i] * h.digestionduration
        if Δt > timeuntildigested
            newΔt = Δt - timeuntildigested
            Δt = timeuntildigested
            h.stomach[i] = 0.
            decreasehunger!(h, Δt, i)
            update!(h, newΔt, ismdpresent, i)
        else
            h.stomach[i] -= Δt/h.digestionduration
            decreasehunger!(h, Δt, i)
        end
    else
        increasehunger!(h, Δt, ismdpresent, i)
    end
end

struct SimpleHunger{C,T}
    hungertimeconstant::T
    digestiontimeconstant::T
    hunger::Vector{Float64}
    cachemodulation::C
end
hunger(h::SimpleHunger) = h.hunger
Base.eachindex(h::SimpleHunger) = eachindex(h.hunger)
function eat!(h::SimpleHunger, nutritionvalue::Number, i)
    h.hunger[i] -= h.hunger[i] * nutritionvalue
end
function update!(h::SimpleHunger, Δt, ismdpresent, i)
    increasehunger!(h, Δt, ismdpresent, i)
    update!(h.cachemodulation, Δt/h.hungertimeconstant, i)
end
function duration_above_threshold(h::SimpleHunger, Δt, θ, ismdpresent, i)
    if ismdpresent
        duration_above_threshold_decr(h.hunger[i], θ, Δt, h.digestiontimeconstant)
    else
        duration_above_threshold_incr(h.hunger[i], θ, Δt, h.digestiontimeconstant)
    end
end


struct ModulatedSpecSatParams{N}
    weights::NTuple{N, Float64}
    bias::Float64
end
struct SpecSatParams{N}
    weights::NTuple{N, Float64}
end
struct SpecSatOrthoParams{E, C, N}
    eatpreference::E
    cachepreference::C
    nutritionvalues::NTuple{N,Float64}
    foodtypes::Vector{Food}
    lookup::NTuple{N_FOODTYPE,Int}
end
lookup(s::SpecSatOrthoParams, item) = s.lookup[foodindex(item)]
NestedStructInitialiser.isleaf(::AbstractArray{<:Food}) = true

function getweight(h, params::ModulatedSpecSatParams, i)
    h[i] * params.weights[i] + params.bias
end
function getweight(::Any, params::SpecSatParams, i)
    params.weights[i]
end
function getcacheweight(::HungerModulatedCaching, h, specsat::SpecSatOrthoParams, typ::Int)
    i = lookup(specsat, typ)
    if typ == Int(Stone)
        specsat.cachepreference.weights[i]
    else
        getweight(hunger(h), specsat.cachepreference, i)
    end
end
function getcacheweight(c::CacheModulatedCaching, ::Any, specsat::SpecSatOrthoParams, typ)
    getweight(c.cachemotivation, specsat.cachepreference, lookup(specsat, typ))
end
function getcacheweight(c::CacheModulatedCaching2, ::Any, specsat::SpecSatOrthoParams, typ)
    getweight(c.hunger, specsat.cachepreference, lookup(specsat, typ))
end
function eat!(h, specsat::SpecSatOrthoParams, item)
    (item.freshness == 0 || item.id == Stone) && return nothing
    i = lookup(specsat, item)
    eat!(h, specsat.nutritionvalues[i], i)
end

abstract type SpecSatAgents end
struct SpecSatAgent{Ts, Th} <: SpecSatAgents
    params::Params
    specsatparams::Ts
    inspectpreference::Float64
    hungermodel::Th
end

function eatfromweight(agent::SpecSatAgents, object)
    geteatweight(agent, foodindex(object))
end
function cachefromtoweight(agent::SpecSatAgent, object, ::Any)
    getcacheweight(agent, foodindex(object))
end
inspectweight(agent::SpecSatAgent, ::Any) = agent.inspectpreference

function geteatweight(agent::SpecSatAgents, typ::Int)
    s = agent.specsatparams
    typ == Int(Stone) && return 0.
    getweight(hunger(agent.hungermodel), s.eatpreference, lookup(s, typ))
end

function eatcallback!(agent::SpecSatAgents, item, trays)
    eat!(agent.hungermodel, agent.specsatparams, item)
	abundancedecrease!(agent, item.id, trays)
end
abundancedecrease!(::Any, ::Any, ::Any) = nothing
function cachecallback!(agent::SpecSatAgents, item, tray)
    typ = item.id
    add!(tray, typ, 1)
    modifycachemotivation!(agent.hungermodel, agent.specsatparams, typ)
end
modifycachemotivation!(::Any, ::Any, ::Any) = nothing


function getcacheweight(agent::SpecSatAgents, typ)
    h = agent.hungermodel
    getcacheweight(h.cachemodulation, h, agent.specsatparams, typ)
end

function updateagentstate!(agent::SpecSatAgent, cage, Δt)
    updatehunger!(agent, cage, Δt)
end
updatehunger!(agent, cage, Δt) = update!(agent.hungermodel, Δt, cage.ismdpresent)

