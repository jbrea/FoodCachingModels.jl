mutable struct Cage{T}
	t::T
    timeunit::T
	ismdpresent::Bool
	inspectionobs::Array{InspectionObserver, 1}
	cacheableitems::Array{FoodItem, 1}
	eatableitems::Array{FoodItem, 1}
	trays::Array{Tray, 1}
    opentrays::Array{Tray, 1}
end
Cage(; timeunit = 1.0u"minute", t = 0.0u"minute", kargs...) =
    Cage(t, timeunit, false,
                      InspectionObserver[], FoodItem[], FoodItem[], Tray[], Tray[])

struct Model{integrationmode, Agent, T, Tracker}
	agent::Agent
    cage::Cage{T}
	tracker::Tracker
end
function Model{I}(agent::A, cage::Cage{T}, tracker::B) where {I, A, T, B}
    Model{I, A, T, B}(agent, cage, tracker)
end

add!(m::Model, ::Type{MaintenanceDiet}) = m.cage.ismdpresent = true
function add!(m::Model, obj::FoodItem)
	obj.eatable && add!(m.cage.eatableitems, obj)
    if obj.cacheable
	    add!(m.cage.cacheableitems, obj, only_push = obj.eatable)
    end
end
add!(m::Model, tray::Tray) = push!(m.cage.trays, tray)
add!(m::Model, obj::InspectionObserver) = push!(m.cage.inspectionobs, obj)
function checkopentrays!(cage::Cage)
    empty!(cage.opentrays)
    for tray in cage.trays
        tray.closed && continue
        push!(cage.opentrays, tray)
    end
end

remove!(m::Model, ::Type{MaintenanceDiet}) = m.cage.ismdpresent = false
function remove!(m::Model, ::Type{Any})
    m.cage.ismdpresent = false
    for listname in (:trays, :opentrays,
                     :cacheableitems, :eatableitems,
                     :inspectionobs)
        empty!(getfield(m.cage, listname))
    end
end
remove!(m::Model, obj) = removefromobjlists!(m, obj)
function remove!(m::Model, ::Tray)
    empty!(m.cage.trays)
    empty!(m.cage.opentrays)
end
function removefromlist!(list::Array, obj)
	ind = findfirst(x -> x == obj, list)
	splice!(list, ind)
end
removefromobjlists!(m::Model, t) = nothing
function removefromobjlists!(m::Model, tray::Tray)
	removefromlist!(m.cage.trays, tray)
    if !tray.closed
	    removefromlist!(m.cage.opentrays, tray)
    end
end
removefromobjlists!(m::Model, obj::FoodItem) = removefromobjlists!(m.cage, obj)
function removefromobjlists!(container::T, obj::FoodItem) where T
    if obj.n == 1
        if T <: Cage && obj.cacheable
            removefromlist!(container.cacheableitems, obj)
        end
        if obj.eatable || T <: Tray
            removefromlist!(container.eatableitems, obj)
        end
    else
        obj.n -= 1
    end
    return nothing
end
removefromobjlists!(m::Model, ::InspectionObserver) = empty!(m.cage.inspectionobs)

function selecttray(::Any, trays, item)
	tray = sample(trays)
	if !tray.closed
		return tray
	end
	opentrays = findall(x -> x.closed == false, trays)
    return sample(trays[opentrays])
end

function eatfrom!(agent, object, i, trays)
    item = object.eatableitems[i]
    eatcallback!(agent, item, trays)
    inspectcallback!(agent, object, item)
    removefromobjlists!(object, item)
    return nothing
end
function cachefromto!(agent, tray::Tray, i, trays, j)
    item = tray.eatableitems[i]
    inspectcallback!(agent, tray, item)
    if j == 0 || trays[j] == tray
        return nothing
    else
        cachefromto!(agent, tray, item, trays[j])
    end
    return nothing
end
function cachefromto!(agent, cage::Cage, i, trays, j)
    cachefromto!(agent, cage, cage.cacheableitems[i], trays[j])
end
function cachefromto!(agent, container, item, tray)
    cachecallback!(agent, item, tray)
    removefromobjlists!(container, item)
    return nothing
end
function _inspect!(m::Model, tray, trays)
    if length(m.cage.inspectionobs) > 0 # inspection counter
		push!(m.cage.inspectionobs[1].trayappearances, tray.appearance)
	end
    if length(tray.eatableitems) > 0
        return _act!(m, true, false, true, tray, tray.eatableitems, trays,
                     otherpossible = false)
    else
        inspectcallback!(m.agent, tray) # pilfer
        return 0.
    end
end

function decide!(m::Model)
    eatpossible = length(m.cage.eatableitems) > 0
    inspectpossible = length(m.cage.trays) > 0
    cachepossible = false
    opentrays = m.cage.opentrays
    if inspectpossible && length(m.cage.cacheableitems) > 0
        cachepossible = length(opentrays) > 0
    end
    _act!(m, eatpossible, inspectpossible, cachepossible,
          m.cage, m.cage.cacheableitems, opentrays)
end
function updatecut(cut, wmax, w, t)
    if w > wmax
        wmax = w
    end
    t > 0 ? cut : wmax, wmax
end
function squash(agent)
#     c = exp(agent.params.u)
    w -> begin
        w = min(1., max(eps(), w))
#         c*w/(1 + w*(c - 1))
    end
end
function _act!(m::Model, eatpossible, inspectpossible, cachepossible,
               object, cacheableitems, opentrays; otherpossible = true)
    lcache = cachepossible * length(cacheableitems)
    leat = eatpossible * length(object.eatableitems)
    linspect = inspectpossible * length(opentrays)
    c3 = lcache/(lcache + otherpossible)
    c2 = linspect/(linspect + lcache + otherpossible)
    c1 = leat/(leat + linspect + lcache + otherpossible)
    cut = 1. # max(5, m.agent.params.otheractionweight)
    wmax = otherpossible ? m.agent.params.otheractionweight : eps()
#     @show c1 c2 c3 leat linspect lcache eatpossible inspectpossible cachepossible
    t = 5 * (leat + linspect + lcache)
    lookup = Dict{Tuple{Symbol,Int,Int}, Float64}()
    transfer = squash(m.agent)
    while true
        t -= 1
        if eatpossible && rand() < c1
            i = rand(1:length(object.eatableitems))
            key = (:eat, i, 0)
            if haskey(lookup, key)
                w = lookup[key]
            else
                w = eatfromweight(m.agent, object.eatableitems[i]) |> transfer
                lookup[key] = w
            end
            cut, wmax = updatecut(cut, wmax, w, t)
            if w > cut * rand()
                trackeat!(m.tracker, w, object.eatableitems[i], m.cage.t)
                eatfrom!(m.agent, object, i, m.cage.trays)
                return timeout_eat(m.agent)
            end
        elseif inspectpossible && rand() < c2
            i = rand(1:length(opentrays))
            key = (:inspect, i, 0)
            if haskey(lookup, key)
                w = lookup[key]
            else
                w = inspectweight(m.agent, opentrays[i]) |> transfer
                lookup[key] = w
            end
            cut, wmax = updatecut(cut, wmax, w, t)
            if w > cut * rand()
                trackinspect!(m.tracker, w, opentrays[i], m.cage.t)
                _inspect!(m, opentrays[i], opentrays)
                return timeout_inspect(m.agent)
            end
        elseif cachepossible && rand() < c3
            i_item = rand(1:length(cacheableitems))
            i_tray = rand(1:length(opentrays))
            key = (:cache, i_item, i_tray)
            if haskey(lookup, key)
                w = lookup[key]
            else
                w = cachefromtoweight(m.agent, cacheableitems[i_item], opentrays[i_tray]) |> transfer
                lookup[key] = w
            end
            cut, wmax = updatecut(cut, wmax, w, t)
            if w > cut * rand()
                trackcache!(m.tracker, w, cacheableitems[i_item], opentrays[i_tray], m.cage.t)
                cachefromto!(m.agent, object, i_item, opentrays, i_tray)
                return timeout_cache(m.agent)
            end
        elseif otherpossible
            w = m.agent.params.otheractionweight |> transfer
            if w > cut * rand()
                trackotheraction!(m.tracker, w, m.cage.t)
                return timeout_other(m.agent)
            end
        end
    end
end


elapsed(m::Model) = m.cage.t - m.cage.t0
function wait!(m::Model, Δt, extracondition = m -> false)
    checkopentrays!(m.cage)
	track!(m.tracker, m, m.cage.t)
	tfinal = m.cage.t + Δt
    t = m.cage.t
	while t < tfinal
        if length(m.cage.eatableitems) == 0 && length(m.cage.trays) == 0 || extracondition(m)
            t = tfinal
        else
            t += decide!(m)
        end
        t = min(t, tfinal)
        integrate!(m, t)
    end
end

integrate!(m::Model, t) = _integrate!(m, t)
function _integrate!(m, t)
    updateagentstate!(m.agent, m.cage, t - m.cage.t)
    m.cage.t = t
    track!(m.tracker, m, m.cage.t)
end
function integrate!(m::Model{:euler}, t)
    s = m.cage.t
    while true
        s += 60u"s"
        t < s && return _integrate!(m, t)
        _integrate!(m, s)
    end
end


function countfooditems(m::Model, typ)
	n = countfooditems(m.cage.eatableitems, typ)
	for tray in m.cage.trays
		n += countfooditems(tray.eatableitems, typ)
	end
	n
end

eatcallback!(::Any, ::Any, ::Any) = nothing
inspectcallback!(::Any, ::Any, item = nothing) = nothing
function cachecallback!(::Any, obj, tray)
	add!(tray, obj.id, 1)
end

selecttraystoinspect(::Any, trays) = sample(trays)
updateagentstate!(::Any, ::Any, ::Any) = nothing
