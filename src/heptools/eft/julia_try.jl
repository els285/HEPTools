
using Plots
using EFTfitter
using BAT
using IntervalSets
using Distributions
using DataStructures
using ValueShapes
using FileIO
using HDF5
using NestedSamplers
using JSON

your_samples_object = "../results_ptz_laurynas_ethantest/BATobjects/samples.jld2"
key_of_operator=:cHt

samples = unshaped(your_samples_object)
idx = asindex(maybe_shaped_samples, key_of_operator)
s = flatview(samples.v)[idx, :]
edges = StatsBase.histrange((s, ), StatsBase._nbins_tuple((s, ), (200,)), false)[1]
hist = fit(Histogram, s, FrequencyWeights(samples.weight), edges, closed = false)
uvbd = EmpiricalDistributions.UvBinnedDist(hist)
marg = MarginalDist((idx,), uvbd, varshape(your_samples_object))


