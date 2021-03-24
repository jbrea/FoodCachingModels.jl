mutable struct MemoryContent
    hunger::Vector{Float64}
    rewarded::Bool
    degraded::Bool
    pilfered::Bool
end
MemoryContent(hunger, experience) = MemoryContent(copy(hunger),
                                                  experience == :rewarded,
                                                  experience == :degraded,
                                                  experience == :pilfered)
function change!(m::MemoryContent, hunger, experience)
    m.hunger .*= .5
    @. m.hunger += .5 * hunger
    experience === nothing && return m
    setfield!(m, experience, true)
    m
end
function memoryweight(m, p, i)
    h = m.hunger[i]
    h * p.hungerscale +
    (2 * (h > p.dontcare) - 1) * m.rewarded * p.rewardedscale -
    m.degraded * p.degradedscale -
    m.pilfered * p.pilferedscale
end
mutable struct PlanningAgentParams{S,T}
    staticparams::S
    hungerscale::Float64
    rewardedscale::Float64
    degradedscale::Float64
    pilferedscale::Float64
    discountfactor::Float64
    dontcare::Float64
    timeref::T
    current::Dict{Int, MemoryContent}
    memory::Vector{Dict{Int, MemoryContent}}
end
hasmemory(::Any) = false
hasmemory(p::PlasticCachingAgent) = hasmemory(p.cacheparams)
hasmemory(::PlanningAgentParams) = true
freshnessweight(p::PlanningAgentParams, h, f) = freshnessweight(p.staticparams, h, f)
inspectbias(s::PlanningAgentParams) = inspectbias(s.staticparams)
inspectslope(s::PlanningAgentParams) = inspectslope(s.staticparams)

function update_model!(agent, tray, experience = nothing)
    # update long term memory
    p = agent.cacheparams
    if agent.t - p.timeref < 1u"hr" # bind together
        if haskey(p.current, tray)
            change!(p.current[tray], hunger(agent), experience)
        else
            p.current[tray] = MemoryContent(hunger(agent), experience)
        end
    else
        push!(p.memory, p.current)
        length(p.memory) > 30 && popfirst!(p.memory)
        p.current = Dict{Int, MemoryContent}()
        p.current[tray] = MemoryContent(hunger(agent), experience)
    end
    p.timeref = agent.t
    nothing
end
function getquery(p)
    l = length(p.memory)
    s = max(1, l - 1)
    query = Vector{Base.KeySet{Int64,Dict{Int64,MemoryContent}}}(undef, l - s + 2)
    j = 1
    for i in max(1, length(p.memory) - 1):length(p.memory)
        query[j] = keys(p.memory[i])
        j += 1
    end
    query[end] = keys(p.current)
    query
end
function querymatching(memory, q)
    res = zeros(length(memory))
    idx = 1
    for j in 1:min(length(memory), length(q))
        for i in 1:j
            res[idx] += q[end-j+i] == keys(memory[i])
        end
        idx += 1
    end
    for j in 1:length(memory)  - length(q)
        for i in eachindex(q)
            res[idx] += q[i] == keys(memory[j+i])
        end
        idx += 1
    end
    res
end
function idxs_bestmatch(memory, q)
    m = querymatching(memory, q)
    mmax = maximum(m)
    mmax > 0 ? findall(==(mmax), m) : Int[]
end
function _cacheweight(agent, cache_feature, item)
    p = agent.cacheparams
    idxs = idxs_bestmatch(p.memory, getquery(p))
    w = 0.
    γ = 1.
    for i in Iterators.take(reverse(idxs), 5)
        similarity = memorysimilarity(p.memory, p.current, i)
        for j in i+1:min(length(p.memory), i+6)
            if haskey(p.memory[j], cache_feature)
                i = agent.specsatparams.lookup[foodindex(item)]
                w += similarity * γ * memoryweight(p.memory[j][cache_feature], p, i)
            end
            γ *= p.discountfactor
        end
    end
    tanh(w/2)
end
function memorysimilarity(m1::MemoryContent, m2::MemoryContent)
    1 - (sum(abs, m1.hunger .- m2.hunger) +
         abs(m1.pilfered - m2.pilfered) +
         abs(m1.degraded - m2.degraded) +
         abs(m1.rewarded - m2.rewarded))/(length(m1.hunger) + 1)
end
function memorysimilarity(memory, current, idx1, idx2 = length(memory) + 1)
    sim = 0.
    for i in 0:min(idx1-1, idx2-1, 2)
        m1 = memory[idx1 - i]
        m2 = idx2 - i > length(memory) ? current : memory[idx2 - i]
        k1 = keys(m1)
        k2 = keys(m2)
        nkeys = 0
        for k in k1
            haskey(m2, k) && (nkeys += 1)
        end
        nkeys == 0 && continue
        tmpsim = 0.
        for k in intersect(k1, k2)
            tmpsim += memorysimilarity(m1[k], m2[k])
        end
        sim += 1/3 * tmpsim/nkeys
    end
    sim
end
function cachefromtoweight(agent::PlasticCachingAgent{<:PlanningAgentParams},
                           item::FoodItem, tray::Tray)
    cache_feature = cache_features(agent.snapshots, tray)
    update_model!(agent, cache_feature)
    witem = getcacheweight(agent, foodindex(item.id))
    p = agent.cacheparams
    _cacheweight(agent, cache_feature, item) + witem
end
function inspectweight(agent::PlasticCachingAgent{<:PlanningAgentParams}, tray::Tray)
    cache_feature = cache_features(agent.snapshots, tray)
    update_model!(agent, cache_feature)
    _inspectweight(agent, tray)
end
