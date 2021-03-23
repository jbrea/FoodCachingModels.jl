module FoodCachingModels

using FoodCachingExperiments, NestedStructInitialiser, Distributions, Unitful,
      DataFrames, RuntimeGeneratedFunctions, Distributed

import FoodCachingExperiments: add!, remove!, wait!, countfooditems
import FoodCachingExperiments: Tray, Food, FoodItem, InspectionObserver,
                               MaintenanceDiet, observe!, N_FOODTYPE, Stone,
                               foodindex, bload, bsave

export Baseline, MotivationalControl, EpisodicLikeMemory, PlasticCaching, ReplayAndPlan,
       save, load

include("trackers.jl")
include("common.jl")
include("baseline.jl")
include("motivational_control.jl")
include("episodic_memory.jl")
include("plastic_caching.jl")
include("replay_and_plan.jl")
include("population.jl")

end
