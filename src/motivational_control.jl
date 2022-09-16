function foodlookup(foodtypes)
    lookup = zeros(Int, N_FOODTYPE)
    for (i, v) in pairs(foodtypes)
        lookup[foodindex(v)] = i
    end
    tuple(lookup...)
end

struct Hunger{T}
    hungertimeconstant::T
    digestiontimeconstant::T
    digestionduration::T
    hunger::Vector{Float64}
    stomach::Vector{Float64}
end

hunger(h::Hunger) = h.hunger
Base.eachindex(h::Hunger) = eachindex(h.hunger)
hunger(agent) = hunger(agent.hungermodel)
eat!(h::Hunger, nutritionvalue::Number, i) = h.stomach[i] += nutritionvalue
function eat!(h, nutritionvalue)
    for i in eachindex(nutritionvalue)
        eat!(h, nutritionvalue[i], i)
    end
end
function duration_above_threshold(h::Hunger, Δt, θ, ismdpresent, i)
    if ismdpresent
        h.hunger[i] <= θ && return 0.0u"s"
        return min(h.digestiontimeconstant * log(h.hunger[i]/θ), Δt)
    end
    if h.stomach[i] > 0.
        timeuntildigested = h.stomach[i] * h.digestionduration
        if Δt > timeuntildigested
            hcopy = deepcopy(h)
            newΔt = Δt - timeuntildigested
            hcopy.stomach[i] = 0.
            d_above = h.hunger[i] <= θ ? 0.0u"s" : min(h.digestiontimeconstant * log(h.hunger[i]/θ), timeuntildigested)
            decreasehunger(hcopy, timeuntildigested, i)
            d_above + duration_above_threshold(hcopy, newΔt, θ, ismdpresent, i)
        else
            h.hunger[i] <= θ && return 0.0u"s"
            return min(h.digestiontimeconstant * log(h.hunger[i]/θ), Δt)
        end
    else
        h.hunger[i] >= θ && return Δt
        return max(0.0u"s", Δt - h.digestiontimeconstant * log(θ/h.hunger[i]))
    end
end
function update!(h, Δt, ismdpresent = false)
    for i in eachindex(h)
        update!(h, Δt, ismdpresent, i)
    end
    hunger(h)
end
decreasehunger(h, Δt, i) = h.hunger[i] *= exp(-Δt/h.digestiontimeconstant)
function increasehunger!(h, Δt, ismdpresent, i)
    ismdpresent && return decreasehunger(h, Δt, i)
    h.hunger[i] = 1 + (h.hunger[i] - 1) * exp(-Δt/h.hungertimeconstant)
end
function update!(h::Hunger, Δt, ismdpresent, i)
    if h.stomach[i] > 0.
        timeuntildigested = h.stomach[i] * h.digestionduration
        if Δt > timeuntildigested
            newΔt = Δt - timeuntildigested
            Δt = timeuntildigested
            h.stomach[i] = 0.
            decreasehunger(h, Δt, i)
            update!(h, newΔt, ismdpresent, i)
        else
            h.stomach[i] -= Δt/h.digestionduration
            decreasehunger(h, Δt, i)
        end
    else
        increasehunger!(h, Δt, ismdpresent, i)
    end
end

struct SimpleHunger{T}
    hungertimeconstant::T
    digestiontimeconstant::T
    hungerdecreasescale::Float64
    hunger::Vector{Float64}
end
hunger(h::SimpleHunger) = h.hunger
Base.eachindex(h::SimpleHunger) = eachindex(h.hunger)
function eat!(h::SimpleHunger, nutritionvalue::Number, i)
    h.hunger[i] -= h.hungerdecreasescale * h.hunger[i] * nutritionvalue
end
update!(h::SimpleHunger, Δt, ismdpresent, i) = increasehunger!(h, Δt, ismdpresent, i)
function duration_above_threshold(h::SimpleHunger, Δt, θ, ismdpresent, i)
    if ismdpresent
        h.hunger[i] ≤ θ && return 0.0u"s"
        return min(h.digestiontimeconstant * log(h.hunger[i]/θ), Δt)
    else
        h.hunger[i] ≥ θ && return Δt
        return max(0.0u"s", Δt - h.digestiontimeconstant * log(θ/h.hunger[i]))
    end
end

struct SimpleFactorizedHunger{T}
    hungertimeconstant::T
    digestiontimeconstant::T
    cachedecreasescale::Float64
    hunger::Vector{Float64}
    cachemotivation::Vector{Float64}
end
hunger(h::SimpleFactorizedHunger) = h.hunger
Base.eachindex(h::SimpleFactorizedHunger) = eachindex(h.hunger)
function eat!(h::SimpleFactorizedHunger, nutritionvalue::Number, i)
    h.hunger[i] -= h.hunger[i] * nutritionvalue
end
function update!(h::SimpleFactorizedHunger, Δt, ismdpresent, i)
    increasehunger!(h, Δt, ismdpresent, i)
    h.cachemotivation[i] = 1 + (h.cachemotivation[i] - 1) * exp(-Δt/h.hungertimeconstant)
end
function modifycachemotivation!(h::SimpleFactorizedHunger, specsat, typ)
    i = lookup(specsat, typ)
    h.cachemotivation[i] -= h.cachemotivation[i] * h.cachedecreasescale
end
function duration_above_threshold(h::SimpleFactorizedHunger, Δt, θ, ismdpresent, i)
    if ismdpresent
        h.hunger[i] ≤ θ && return 0.0u"s"
        return min(h.digestiontimeconstant * log(h.hunger[i]/θ), Δt)
    else
        h.hunger[i] ≥ θ && return Δt
        return max(0.0u"s", Δt - h.digestiontimeconstant * log(θ/h.hunger[i]))
    end
end

struct FactoredMotivationalControl{H,C}
    hungermodulation::H
    cachemodulation::C
end
hunger(h::FactoredMotivationalControl) = hunger(h.hungermodulation)
Base.eachindex(h::FactoredMotivationalControl) = eachindex(h.hungermodulation)
eat!(h::FactoredMotivationalControl, nutritionvalue::Number, i) =
    eat!(h.hungermodulation, nutritionvalue, i)
function update!(h::FactoredMotivationalControl, Δt, ismdpresent, i)
    update!(h.hungermodulation, Δt, ismdpresent, i)
    update!(h.cachemodulation, Δt, false, i)
end
duration_above_threshold(h::FactoredMotivationalControl, Δt, θ, ismdpresent, i) =
    duration_above_threshold(h.hungermodulation, Δt, θ, ismdpresent, i)
function modifycachemotivation!(h::FactoredMotivationalControl, specsat, typ)
    eat!(h.cachemodulation, 1., lookup(specsat, typ))
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
function getcacheweight(h::Union{SimpleHunger, Hunger}, specsat::SpecSatOrthoParams, typ::Int)
    i = lookup(specsat, typ)
    if typ == Int(Stone)
        specsat.cachepreference.weights[i]
    else
        getweight(hunger(h), specsat.cachepreference, i)
    end
end
function getcacheweight(h::FactoredMotivationalControl, specsat::SpecSatOrthoParams, typ)
    getweight(hunger(h.cachemodulation), specsat.cachepreference, lookup(specsat, typ))
end
function getcacheweight(h::SimpleFactorizedHunger, specsat::SpecSatOrthoParams, typ)
    getweight(h.cachemotivation, specsat.cachepreference, lookup(specsat, typ))
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
    getcacheweight(agent.hungermodel, agent.specsatparams, typ)
end

function updateagentstate!(agent::SpecSatAgent, cage, Δt)
    updatehunger!(agent, cage, Δt)
end
updatehunger!(agent, cage, Δt) = update!(agent.hungermodel, Δt, cage.ismdpresent)

