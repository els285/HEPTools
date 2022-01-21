using BAT
using FileIO
using ArraysOfArrays
using HDF5
using DelimitedFiles
using Tables
using TypedTables
using ValueShapes

function return_1d_sample_points(file_name::String,i::Int64)

    full_object = load(file_name)

    samples = full_object["samples"]

    wc_array = Tables.columns(samples.v)[i]

    return wc_array

end


