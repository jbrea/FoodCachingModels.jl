struct Init{T}
    i::T
end
(i::Init)() = deepcopy(i.i)
function defaults(::Type{Cage}; kwargs...)
    merge((t = .0u"minute",
           timeunit = 1.0u"minute",
           ismdpresent = false,
           inspectionobs = Init(InspectionObserver[]),
           cacheableitems = Init(FoodItem[]),
           eatableitems = Init(FoodItem[]),
           trays = Init(Tray[]),
           opentrays = Init(Tray[]),
           tracker = NoTracker),
          kwargs)
end
function defaults(::Type{Model}; kwargs...)
    merge((cage = Cage, agent = PlasticCachingAgent),
          kwargs)
end
function defaults(::Type{SpecSatAgent}; kwargs...)
    merge((hungermodel = Hunger, specsatparams = SpecSatOrthoParams),
          kwargs)
end
default_foodtypes(::Any) = collect(instances(Food))
default_foodtypes(es::Vector{Symbol}) = union(vcat(FoodCachingExperiments.foodtypes.(es)...))
function defaults(::Type{SpecSatOrthoParams};
                  experiments = nothing,
                  foodtypes = default_foodtypes(experiments),
                  cachepreferencemodel = ModulatedSpecSatParams,
                  eatpreferencemodel = ModulatedSpecSatParams,
                  kwargs...)
    n = length(foodtypes)
    n1 = n - (Stone in foodtypes)
    merge(merge((lookup = foodlookup(foodtypes),
                 foodtypes = foodtypes),
                kwargs),
          (specsatparams = SpecSatOrthoParams{eatpreferencemodel{n1}, cachepreferencemodel{n}, n1},))
end
function defaults(::Type{Hunger};
                  timeunit = 1.0u"minute",
                  experiments = nothing,
                  foodtypes = default_foodtypes(experiments),
                  cachemodulation = HungerModulatedCaching,
                  kwargs...)
    n = length(foodtypes)
    merge(merge((hunger = Init(zeros(n)),
                 stomach = Init(zeros(n)),
                 timeunit = timeunit,
                 cachemodulation = cachemodulation,
                 cachemotivation = Init(ones(n)),
                 foodtypes = foodtypes),
                kwargs),
          (hungermodel = Hunger{cachemodulation,typeof(timeunit)},))
end
function defaults(::Type{CacheModulatedCaching2};
                  timeunit = 1.0u"minute",
                  experiments = nothing,
                  foodtypes = default_foodtypes(experiments),
                  kwargs...)
    merge(merge((; timeunit, experiments, foodtypes), kwargs),
          (cachemodulation = CacheModulatedCaching2{typeof(timeunit), length(foodtypes)},))
end
function defaults(::Type{SimpleHunger};
                  timeunit = 1.0u"minute",
                  experiments = nothing,
                  foodtypes = default_foodtypes(experiments),
                  cachemodulation = HungerModulatedCaching,
                  kwargs...)
    n = length(foodtypes)
    merge(merge((hunger = Init(zeros(n)),
                 timeunit = timeunit,
                 cachemodulation = cachemodulation,
                 cachemotivation = Init(zeros(n)),
                 foodtypes = foodtypes),
                kwargs),
          (hungermodel = SimpleHunger{cachemodulation,typeof(timeunit)},))
end
function defaults(::Type{PlasticCachingAgent};
                  kwargs...)
    merge((snapshots = SnapshotSimple,
           fillleveldecreasemodel = NoFillLevelDecrease,
           staticparams = StaticCachingAgentParams{AdditiveFreshness},
           freshnessestimator = PiecewiseConstantFreshnessEstimator,
           cacheparams = PlasticCachingAgentParams,
           compensatorycaching = HungerIncrease,
           trayweights = Init(Dict{Tuple{Int, Food}, Float64}()),
          ),
          defaults(SpecSatAgent; kwargs...))
end
function defaults(::Type{SnapshotSimple};
                  timeunit = 1.0u"minute",
                  snapshotfeatures = AppearanceFeature,
                  kwargs...)
    K = snapshotfeatures == BothFeatures ? Tuple{Int, Int} : Int
    T = typeof(timeunit)
    C = SnapshotContent{T}
    merge(merge((mapping = Init(Dict{Int, Dict{Food, C}}()),
                 timeunit = timeunit,
                 snapshotfeatures = snapshotfeatures),
                kwargs),
          (snapshots = SnapshotSimple{snapshotfeatures, K, T},))
end
function defaults(::Type{HungerIncrease}; timeunit = 1.0u"minute", kwargs...)
    merge(merge((timeunit = timeunit,), kwargs),
          (compensatorycaching = HungerIncrease{typeof(timeunit)},))
end
function defaults(::Type{PiecewiseConstantFreshnessEstimator}; kwargs...)
    merge((freshness = Init([PiecewiseConstant([1u"d", 2u"d", 3u"d", 4u"d", 7u"d", 15.0u"d", Inf*u"d"], [1., 1., 1., 1., 1., .5, 0.]) for _ in 1:N_FOODTYPE]),),
          kwargs)
end
function defaults(::Type{PlanningAgentParams}; kwargs...)
    merge((timeref = .0u"minute", current = Init(Dict{Int, MemoryContent}()),
           memory = Init(Dict{Int, MemoryContent}[])), kwargs)
end
defaults(::Any; kwargs...) = kwargs
function defaults(; model = Model, kwargs...)
    def = merge((model = model,), kwargs)
    i = 1
    while i <= length(def)
        def = merge(def, defaults(def[i]; def...))
        i += 1
    end
    def
end

const BOUNDS = Dict{Symbol, Tuple{Float64, Float64}}()
# SimpleAgent
BOUNDS[:otheractionweight] = (0., 1.)
BOUNDS[:timeout_eat] = (10., 200.)
BOUNDS[:timeout_cache] = (10., 200.)
BOUNDS[:timeout_inspect] = (10., 200.)
BOUNDS[:timeout_other] = (10., 200.)
BOUNDS[:eatpreference] = (0., 4.)
BOUNDS[:cachepreference] = (0., 4.)
BOUNDS[:inspectpreference] = (0., 4.)
# SpecSatAgent
BOUNDS[:hungertimeconstant] = (50., 300.)
BOUNDS[:digestiontimeconstant] = (eps(), 20.)
BOUNDS[:digestionduration] = (.5, 10.)
BOUNDS[:weights] = (0., 1.)
BOUNDS[:bias] = (-1., 1.)
BOUNDS[:nutritionvalues] = (.1, 1.)
BOUNDS[:update_value] = (.1, 1.)
BOUNDS[:inspectpreference] = (0., 4.)
# PlasticCachingAgent
BOUNDS[:initialweight] = (0., 1.)
BOUNDS[:inspectbias] = (-2., .5)
BOUNDS[:inspectslope] = (0., 5.)
BOUNDS[:hungerincreasethreshold] = (.9, 1.)
BOUNDS[:hungerincreasetimeconstant] = (100., 500.)
BOUNDS[:freshness_alpha] = (0., 1.)
BOUNDS[:maxcacheweight] = (.2, 1.)
BOUNDS[:pilferdecrease] = (.8, 1.)
BOUNDS[:degradedecrease] = (.8, 1.)
BOUNDS[:retrievalrewardfactor] = (.1, 1.)
BOUNDS[:dontcare] = (0., .5)
# ReplayAndPlan
BOUNDS[:hungerscale] = (0., 2.)
BOUNDS[:rewardedscale] = (0., 1.)
BOUNDS[:degradedscale] = (0., 1.)
BOUNDS[:pilferedscale] = (0., 1.)
BOUNDS[:discountfactor] = (0., 1.)
BOUNDS[:cachedecreasescale] = (0., 1.)

function bounds(p::NestedStructInitialiser.Parameters; kwargs...)
    b = merge(BOUNDS, kwargs)
    l = Float64[]
    u = Float64[]
    for (x, t) in p.free
        bx = b[x]
        for _ in 1:NestedStructInitialiser.free_param_length(t)
            push!(l, bx[1])
            push!(u, bx[2])
        end
    end
    l, u
end

mutable struct Population{T,D}
    m::Vector{Float64}
    s::T
    l::Vector{Float64}
    u::Vector{Float64}
    p::NestedStructInitialiser.Parameters
    dist::D
    constructor
end
function Base.show(io::IO, ::MIME"text/plain", pop::Population)
    println(io, "Population of $(pop.p.fixed[1][2]) ($(pop.dist))")
end
function init(p; pids = procs())
    @distributed vcat for _ in pids
        [initialiser(p)]
    end
end
function setparameters!(m::Population{<:AbstractVector}, x)
    n = div(length(x), 2)
    m.m .= @view x[1:n]
    m.s .= @view x[n+1:end]
    m
end
function setparameters!(m::Population{<:Nothing}, x)
    m.m .= x
    m
end
function Population(; integrationmode = :eventbased,
                      distribution = truncnorm,
                      pids = procs(),
                      kwargs...)
    p = parameters(Model{integrationmode};
                   defaults(; kwargs...)...)
    l, u = bounds(p)
    c = init(p; pids)
    Population(random_init(distribution, length(l))..., l, u, p, distribution, c)
end
softplus(x) = x > 40 ? x : log(exp(x) + 1) + eps()
delta() = nothing
beta(d, s) = d < 0 ? Beta(softplus(s), softplus(s + d)) :
                     Beta(softplus(s - d), softplus(s))
truncnorm(m, s) = TruncatedNormal(m, softplus(s), 0, 1)
random_init(::typeof(beta), l) = (m = 20*rand(l) .- 10, s = 50 * rand(l) .+ 75)
random_init(::typeof(truncnorm), l) = (m = rand(l), s = 4*rand(l) .- 4)
random_init(::typeof(delta), l) = randn(l), nothing
function Base.rand(m::Population{<:AbstractVector})
    m.constructor[myid()](@. (m.u - m.l) * rand(m.dist(m.m, m.s)) + m.l)
end
function Base.rand(m::Population{<:Nothing})
    m.constructor[myid()](@. (m.u - m.l) * 1/(1 + exp(m.m)) + m.l)
end
Base.rand(m::Population, n::Int) = [rand(m) for _ in 1:n]

function save(path, m::Population, dict = Dict{Symbol, Any}())
    dict[:model] = m
    bsave(path, dict)
end
function load(path; pids = procs())
    m = bload(path, @__MODULE__)[:model]
    Population(m.m, m.s, m.l, m.u, m.p, m.dist, init(m.p; pids))
end

Baseline(; kwargs...) = Population(; agent = SimpleAgent, kwargs...)
MotivationalControl(; kwargs...) = Population(; agent = SpecSatAgent,
                                              digestionduration = 5.0u"minute",
                                              kwargs...)
EpisodicLikeMemory(; kwargs...) = Population(; maxcacheweight = 1.,
                                             digestionduration = 5.0u"minute",
                                             cacheparams = StaticCachingAgentParams{AdditiveFreshness},
                                             initialweight = 0., kwargs...)
PlasticCaching(; kwargs...) = Population(; digestionduration = 5.0u"minute",
                                         kwargs...)
ReplayAndPlan(; kwargs...) = Population(; cacheparams = PlanningAgentParams,
                                        digestionduration = 5.0u"minute",
                                        kwargs...)
