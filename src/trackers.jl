struct Tracker
	fieldnames::Array{Any, 1}
	data::DataFrame
end
struct NoTracker end

function Tracker(; trackedfields = [], kargs...)
	if length(trackedfields) == 0
		NoTracker()
	else
		Tracker(trackedfields, DataFrame(t = Quantity[], field = [], value = []))
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

track!(::NoTracker, w, t) = nothing
function track!(tracker::Tracker, m, t)
	for field in tracker.fieldnames
        field ∉ (:weights, :actions) && push!(tracker.data, [t, field, getdeepfield(m, field)])
	end
end

trackweights!(::NoTracker, w, t) = nothing
function trackweights!(tracker::Tracker, w, t)
	if :weights in tracker.fieldnames
        push!(tracker.data, [t, :weights, deepcopy(w)])
	end
end

trackeat!(::NoTracker, w, o, t) = nothing
function trackeat!(tracker::Tracker, w, o, t)
    if :actions in tracker.fieldnames
        push!(tracker.data, [t, :actions, ("eat", w, o.id)])
    end
end
trackinspect!(::NoTracker, w, o, t) = nothing
function trackinspect!(tracker::Tracker, w, tray, t)
    if :actions in tracker.fieldnames
        push!(tracker.data, [t, :actions, ("inspect", w, tray.appearance)])
    end
end
trackcache!(::NoTracker, w, o, tray, t) = nothing
function trackcache!(tracker::Tracker, w, o, tray, t)
    if :actions in tracker.fieldnames
        push!(tracker.data, [t, :actions, ("cache", w, o.id, tray.appearance)])
    end
end
trackotheraction!(::NoTracker, w, t) = nothing
function trackotheraction!(tracker::Tracker, w, t)
    if :actions in tracker.fieldnames
        push!(tracker.data, [t, :actions, ("other", w)])
    end
end

