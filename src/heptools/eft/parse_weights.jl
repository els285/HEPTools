
function return_weights(file_name::String)

    full_object = load(file_name)

    samples = full_object["samples"]

    weights = samples.weight

    return weights

end
