struct Params
    otheractionweight::Float64
    timeout_eat::Float64
    timeout_cache::Float64
    timeout_inspect::Float64
    timeout_other::Float64
end
struct SimpleAgentParams
	eatpreference::Float64
    cachepreference::Float64
	inspectpreference::Float64
end
struct SimpleAgent
    params::Params
	simpleparams::SimpleAgentParams
end

eatfromweight(agent::SimpleAgent, ::Any) = agent.simpleparams.eatpreference
cachefromtoweight(agent::SimpleAgent, ::Any, ::Any) = agent.simpleparams.cachepreference
inspectweight(agent::SimpleAgent, ::Any) = agent.simpleparams.inspectpreference

timeout(T) = (((T - 1) * rand() + 1)/60)u"minute"
timeout_eat(agent) = timeout(agent.params.timeout_eat)
timeout_cache(agent) = timeout(agent.params.timeout_cache)
timeout_inspect(agent) = timeout(agent.params.timeout_inspect)
timeout_other(agent) = timeout(agent.params.timeout_other)
