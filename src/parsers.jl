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
    cut_mask = Vector{String}(undef, file_vector.size[1])
    nml_name = String
    for (i, line) in enumerate(file_vector)
        if match(r"^&", line) !== nothing
            # Double split for whitespace-safety
            nml_name = split(split(line, "&")[2])[1]
            cut_mask[i] = "head " * nml_name
        else
            cut_mask[i] = nml_name
        end
    end

    # First lets go about the allocations as the constructor depends on it
    numelements = 0
    numatoms = 0
    for line in file_vector[cut_mask .== "allocations"]
        nml_attr = split(split(line, "=")[1])[1]
        nml_val = parse(Int64, split(split(line, "=")[2])[1])
        if nml_attr == "numelements"
            numelements = nml_val
        elseif nml_attr == "numatoms"
            numatoms = nml_val
        else
            Logging.@error "Attribute defined in allocations that is not valid!"
        end
    end

    # Finally initialize the empty ebInputData Object
    dataobject = ebInputData(input_file_path, numelements, numatoms)

    # Continue with the crystal info
    for (i, line) in enumerate(file_vector[cut_mask .== "crystal_info"])
        # Disect the line first into attribute, (index, )value
        disect = disect_nml_line(line)

        # Get the type from the attribute out of disect
        attr = disect[1]
        ele_type, conc_type = type_from_attr(attr, dataobject.crystal_info)

        if isa(dataobject.crystal_info[attr], Array)
            # val = split(disect[2])
            println("ARRAY ", attr, " ", conc_type, " ", ele_type, " ", disect[2])
        else
            clean_val = clean_simple_val(disect[2], conc_type)
            println("SIMPLE ", attr, " ", conc_type, " ", ele_type, " ", clean_val)
            dataobject.crystal_info[attr] = parse(conc_type, clean_val)
            Logging.@debug "Updated crystal_info" dataobject.crystal_info
        end
    end
end

function clean_simple_val(val::AbstractString, valtype::DataType)
    if valtype == String || valtype == Bool
        clean_val = strip(val, ['"', ''', '.', ' '])
        return clean_val
    else
        # Other types dont need to be cleaned
        return val
    end
end

function type_from_attr(attr::AbstractString, nml_dict::AbstractDict)
    conc_type = typeof(nml_dict[attr])
    ele_type = eltype(nml_dict[attr])
    return ele_type, conc_type
end

function disect_nml_line(line::AbstractString)
    # First split away the comment
    noco = split(line, "!")[1]
    Logging.@debug "Line without comments" noco

    # Split name-side and value-side
    name = split(split(noco, "=")[1])[1]
    val = split(noco, "=")[2]
    Logging.@debug "Name and Value sides" name val

    # Split the name into attribute and index
    attr_and_index = split(name, "(")
    attr = attr_and_index[1]
    Logging.@debug "Attribute of the Name" attr
    if attr_and_index.size[1] > 1
        index = chop(split(name, "(")[2], tail = 1)
        Logging.@debug "chopped Index of the Name" index

        #Split index into indices
        indices = split(index, ",")
        Logging.@debug "Indices/Ranges of the Name" indices
        return [attr, val, indices]
    else
        return [attr, val]
    end
end

function init_allocations!(
    allocations::AbstractDict,
    numelements::Integer,
    numatoms::Integer,
)
    allocations["numelements"] = numelements
    allocations["numatoms"] = numatoms
end

function init_crystal_info!(ci::AbstractDict, numelements::Integer, numatoms::Integer)
    ci["name"] = "Crystal"
    ci["elements"] = Vector{String}(undef, numelements)
    ci["atomtypes"] = Vector{Int64}(undef, numatoms)
    ci["masses"] = Vector{Float64}(undef, numelements)
    ci["VCA"] = false
    ci["DIB"] = false
    ci["lattvecs"] = Array{Float64}(undef, (3, 3))
    ci["basis"] = Array{Float64}(undef, (3, numatoms))
    ci["polar"] = false
    ci["born"] = Array{Float64}(undef, (3, 3, numatoms))
    ci["epsilon"] = Array{Float64}(undef, (3, 3))

    # This one is falsely labeled on the elphbolt github page
    ci["read_epsiloninf"] = false

    ci["epsiloninf"] = 0.0
    ci["epsilon0"] = 0.0
    ci["T"] = -1.0
    ci["twod"] = false
    ci["subs_masses"] = Vector{Float64}
    ci["subs_conc"] = Vector{Float64}(undef, numelements)
    ci["bound_length"] = 1e12
    ci["specfac"] = 0.0
end

# End of module Parsers
end
