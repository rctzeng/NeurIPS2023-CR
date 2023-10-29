## Best Arm Identification with Fixed Budget: A Large Deviation Perspective
This is the official repository for "Best Arm Identification with Fixed Budget: A Large Deviation Perspective", NeurIPS 2023

### Package packages
 [StatsBase](https://github.com/JuliaStats/StatsBase.jl), 
 [JLD2.jl](https://github.com/JuliaIO/JLD2.jl), 
 [IterTools](https://github.com/JuliaCollections/IterTools.jl), 
 [StatsPlots.jl](https://github.com/JuliaPlots/StatsPlots.jl), 
 [CSV](https://github.com/JuliaData/CSV.jl), 
 [DataFrames](https://github.com/JuliaData/DataFrames.jl)


### Experiments
|Setting |Description                 |
|:------:|:--------------------------:|
|1       |one group of suboptimal arms|
|2       |two group of suboptimal arms|
|3       |linear                      |
|4       |concave                     |
|5       |convex                      |
|6       |stair                       |

We fix `K∈{10,20,40}` for `setting_id∈{1,...,5}` and `K∈{15,21,55}` for `setting_id=6`, while varying the budget `T`.
Each experiment is repeated `N=40,000` times and the averaged error probability `δ` is reported.

### Instructions
See the comments in [experiments.jl](experiments.jl) and [viz.jl](viz.jl) for more details.
 * Run the experiment: `julia -O3 -p8 experiment.jl {setting_id} {M}`, where `-p8` means to use 8 processors
 * Plot the results: `julia -O3 viz.jl {setting_id} {M}`

### Baselines
|Name                                             |Abbrev.  |Description                              |
|:-----------------------------------------------:|:-------:|:---------------------------------------:|
|Continuous Rejects with aggressive rate          |CR-A     |                                         |
|Continuous Rejects with conservative rate        |CR-C     |                                         |
|Successive Rejects                               |SR       |(Audibert, et al., COLT'10)              |
|Sequential Halving                               |SH       |(Karnin, et al., ICML'13)                |
|Uniform Sampling                                 |RR       |                                         |