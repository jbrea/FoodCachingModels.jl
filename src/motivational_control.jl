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
hunger(agent) = hunger(agent.hungermodel)
eat!(h::Hunger, nutritionvalue, i) = h.stomach[i] += nutritionvalue
function eat!(h::Hunger, nutritionvalue)
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
function update!(h::Hunger, Δt, ismdpresent = false)
    for i in 1:length(h.stomach)
        update!(h, Δt, ismdpresent, i)
    end
    h.hunger
end
decreasehunger(h, Δt, i) = h.hunger[i] *= exp(-Δt/h.digestiontimeconstant)
function increasehunger(h, Δt, ismdpresent, i)
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
        increasehunger(h, Δt, ismdpresent, i)
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
NestedStructInitialiser.isleaf(::AbstractArray{<:Food}) = true

function getweight(h::Hunger, params::ModulatedSpecSatParams, i)
    h.hunger[i] * params.weights[i] + params.bias
end
function getweight(::Any, params::SpecSatParams, i)
    params.weights[i]
end
function getcacheweight(h::Hunger, specsat::SpecSatOrthoParams, typ)
    if typ == Int(Stone)
        specsat.cachepreference.weights[specsat.lookup[typ]]
    else
        getweight(h, specsat.cachepreference, specsat.lookup[typ])
    end
end
function eat!(h::Hunger, specsat::SpecSatOrthoParams, item)
    (item.freshness == 0 || item.id == Stone) && return nothing
    i = specsat.lookup[foodindex(item)]
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

function geteatweight(agent::SpecSatAgents, typ)
    s = agent.specsatparams
    typ == Int(Stone) && return 0.
    getweight(agent.hungermodel, s.eatpreference, s.lookup[typ])
end

function eatcallback!(agent::SpecSatAgents, item, trays)
    i = foodindex(item)
    eat!(agent.hungermodel, agent.specsatparams, item)
	abundancedecrease!(agent, item.id, trays)
end
abundancedecrease!(agent, typ, trays) = nothing

function getcacheweight(agent::SpecSatAgents, typ)
    getcacheweight(agent.hungermodel, agent.specsatparams, typ)
end

function updateagentstate!(agent::SpecSatAgent, cage, Δt)
    updatehunger!(agent, cage, Δt)
end
updatehunger!(agent, cage, Δt) = update!(agent.hungermodel, Δt, cage.ismdpresent)

