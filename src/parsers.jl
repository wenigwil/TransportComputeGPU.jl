module Parsers

import Logging

struct ebInputData
    input_file::AbstractString
    allocations::AbstractDict
    crystal_info::AbstractDict

    function ebInputData(
        # Input file that will be read in
        input_file::AbstractString,
        # elphbolt namelists contain quantities that need to be initialized 
        # with dimensions that depend on the numelements and numatoms
        numelements::Integer,
        numatoms::Integer,
        # TODO: For know this struct implements the allocations and 
        # crystal_info namelists from elphbolt. crystal_info depends in 
        # at least 6 quantities on numelements or numatoms.
        # The electrons namelist depends on numconc and numT! Implement 
        # those
    )
        allocations = Dict{String,Any}()
        crystal_info = Dict{String,Any}()

        # Initialize the allocations namelist with the numelements and 
        # numatoms that have to be passed in any way
        init_allocations!(allocations, numelements, numatoms)

        # Initialize the crystal_info namelist. Quantities in this namelist 
        # depend strictly on numelements and numatoms
        init_crystal_info!(crystal_info, numelements, numatoms)

        # Construct the ebInputData-Object
        new(input_file, allocations, crystal_info)
    end
end

function ebInputData(input_file_path::AbstractString)
    # Check that the input_file actually exists
    if isfile(input_file_path) == false
        Logging.@error "Given path is not a regular file!"
        return
    else
        Logging.@debug "Given path is a regular file" input_file_path
    end

    # First save the whole file as a vector of strings
    file_vector = readlines(input_file_path)

    # Cleaning of blank lines and lines of beginning /
    clean_mask = zeros(Int16, file_vector.size[1])
    for (i, line) in enumerate(file_vector)
        if match(r"^/", line) !== nothing
            clean_mask[i] = true
        elseif match(r"^$", line) !== nothing
            clean_mask[i] = true
        end
    end
    # Invert the mask
    clean_mask = convert.(Bool, abs.(clean_mask .- 1))
    # Apply mask
    file_vector = file_vector[clean_mask]

    # Cut at beginning &
    cut_mask = zeros(Int16, file_vector.size[1])
    # TODO implement the array cutting
end

function cut(array_in::AbstractArray, mask::BitArray) end

function init_crystal_info!(
    ci::AbstractDict,
    numelements::Integer,
    numatoms::Integer,
)
    ci["name"] = String
    ci["elements"] = Vector{String}(undef, numelements)
    ci["atomtypes"] = Vector{Int64}(undef, numatoms)
    ci["masses"] = Vector{Float64}(undef, numelements)
    ci["VCA"] = Bool
    ci["DIB"] = Bool
    ci["lattvecs"] = Array{Float64}(undef, (3, 3))
    ci["basis"] = Array{Float64}(undef, (3, numatoms))
    ci["polar"] = Bool
    ci["born"] = Array{Float64}(undef, (3, 3, numatoms))
    ci["epsilon"] = Array{Float64}(undef, (3, 3))

    # This one is falsely labeled on the elphbolt github page
    ci["read_epsiloninf"] = Bool

    ci["epsiloninf"] = Float64
    ci["epsilon0"] = Float64
    ci["T"] = Float64
    ci["twod"] = Bool
    ci["subs_masses"] = Vector{Float64}
    ci["subs_conc"] = Vector{Float64}(undef, numelements)
    ci["bound_length"] = Float64
    ci["specfac"] = Float64
end

function init_allocations!(
    allocations::AbstractDict,
    numelements::Integer,
    numatoms::Integer,
)
    allocations["numelements"] = numelements
    allocations["numatoms"] = numatoms
end

# End of module Parsers
end
