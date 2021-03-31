struct Tracker
	fieldnames::Array{Any, 1}
	data::Dict
    last::Dict
end
struct NoTracker end
struct TimeSeries{T,V}
    t::T
    v::V
end

function Tracker(; trackedfields = [], kargs...)
	if length(trackedfields) == 0
		NoTracker()
	else
        Tracker(trackedfields, Dict(), Dict())
	end
end

function getdeepfield(o, f::Tuple)
	v = getfield(o, f[1])
	if length(f) > 1
		for sf in f[2:end]
			v = getfield(v, sf)
		end
	end
	deepcopy(v)
end
getdeepfield(o, f::Function) = f(o)

function dpush!(dict, key, t, v)
    if haskey(dict, key)
        push!(dict[key].t, t)
        push!(dict[key].v, v)
    else
        dict[key] = TimeSeries([t], [v])
    end
end

track!(::NoTracker, w, t) = nothing
function track!(tracker::Tracker, m, t)
	for field in tracker.fieldnames
        if field ∉ (:actions, )
            v = getdeepfield(m, field)
            k = last(field)
            if haskey(tracker.last, k)
                if !haskey(tracker.data, k) ||
                   !(v == tracker.last[k].v == last(tracker.data[k].v))
                    dpush!(tracker.data, k, tracker.last[k].t, tracker.last[k].v)
                end
            end
            tracker.last[k] = (t = t, v = v)
        end
	end
end

trackeat!(::NoTracker, w, o, t) = nothing
function trackeat!(tracker::Tracker, w, o, t)
    if :actions in tracker.fieldnames
        dpush!(tracker.data, :actions_eat, t, (w, o.id))
    end
end
trackinspect!(::NoTracker, w, o, t) = nothing
function trackinspect!(tracker::Tracker, w, tray, t)
    if :actions in tracker.fieldnames
        dpush!(tracker.data, :actions_inspect, t, (w, tray.appearance))
    end
end
trackcache!(::NoTracker, w, o, tray, t) = nothing
function trackcache!(tracker::Tracker, w, o, tray, t)
    if :actions in tracker.fieldnames
        dpush!(tracker.data, :actions_cache, t, (w, o.id, tray.appearance))
    end
end
trackotheraction!(::NoTracker, w, t) = nothing
function trackotheraction!(tracker::Tracker, w, t)
    if :actions in tracker.fieldnames
        dpush!(tracker.data, :actions_other, t, (w,))
    end
end

