module Phonon

import Logging

"""
    read_ifc2(ifc2_path::String, ebinput_path::String)

Read a Quantum-Espresso-ifc2-output file into a force constant tensor.
"""
function read_qe_ifc2(ifc2_path::String)
    # Check that file exists
    if isfile(ifc2_path) == false
        Logging.@error "Given path is not a regular file!"
        return
    else
        Logging.@debug "Given path is a regular file" ifc2_path
    end

    # Reading the header
    ifc2_file = open(ifc2_path)

    # FIRST LINE
    # Reading number of species (ntype), number of atoms (nat), type of 
    # bravais lattice (?) (ibrav), and a vector of cell dimensions (celldm)
    celldm = Vector{Float64}(undef, 6)
    ntype,
    nat,
    ibrav,
    celldm[1],
    celldm[2],
    celldm[3],
    celldm[4],
    celldm[5],
    celldm[6] = parse.(Float64, split(readline(ifc2_file)))
    # Convert to Ints so we can loop over them
    ntype = convert(Int16, ntype)
    nat = convert(Int16, nat)
    ibrav = convert(Int16, ibrav)

    Logging.@debug "readifc2() First line read..." nat ntype ibrav celldm

    # TODO Deal with the special ibrav=0 case
    # For this we have to find out what quantum espresso prints into the 
    # ifc2-output file in that case... 15 year old forums may suggest a 
    # solution lol

    # NEXT ntype LINE(S)
    # ntype are the number species of atoms
    mass = Vector{Float64}(undef, ntype)
    label = Array{String}(undef, ntype)
    for i = 1:ntype
        temp = split(readline(ifc2_file))
        label[i] = chop(temp[2], head = 1, tail = 1)
        mass[i] = parse(Float64, temp[3])
    end

    Logging.@debug "readifc2() Next $ntype line(s) read..." mass label

    # NEXT nat LINE(S)
    # From these lines we get which species (tipo) sits on which basis 
    # position (basis_pos)
    tipo = Vector{Int16}(undef, nat)
    basis_pos = Matrix{Float64}(undef, (nat, 3))
    for i = 1:nat
        temp = split(readline(ifc2_file))
        tipo[i] = parse(Int16, temp[2])
        basis_pos[i, :] = parse.(Float64, temp[3:end])
    end

    Logging.@debug "readifc2() Next $nat line(s) read..." tipo basis_pos

    # NEXT LINE
    # Is the material a polar one?
    polar_key = split(readline(ifc2_file))[1]

    Logging.@debug "readifc2() Next 1 line read..." polar_key

    if polar_key == "T"
        # NEXT 3 LINES
        # Read the dielectric tensor
        dielectric_tensor = Matrix{Float64}(undef, (3, 3))
        for i = 1:3
            dielectric_tensor[i, :] =
                parse.(Float64, split(readline(ifc2_file)))
        end
        Logging.@debug "readifc2() Next 3 lines read... Dielectric Tensor" dielectric_tensor

        # NEXT 4*nat LINES
        # Read the Born effective charges
        zeff = Array{Float64}(undef, (nat, 3, 3))
        for i = 1:nat
            # We are not reading the first index. We already know that
            readline(ifc2_file)
            for j = 1:3
                zeff[i, j, :] = parse.(Float64, split(readline(ifc2_file)))
            end
        end
        Logging.@debug "readifc2() Next $(4*nat) lines read... Born-Eff. Charges" zeff
    end

    # NEXT LINE
    # Read the super cell sizes
    qcell = parse.(Int32, split(readline(ifc2_file)))
    Logging.@debug "readifc2() Next 1 line read... Super Cell Size" qcell

    # NEXT 3*3*nat*nat*( qcell[1]*qcell[2]*qcell[3] + 1 ) lines
    # First we define the number of force-constant-supercell-fields 
    nfc2 = 3 * 3 * nat * nat
    no_qcells = qcell[1] * qcell[2] * qcell[3]

    # Allocating the force tensor
    ifc2_tensor_dims = (3, 3, nat, nat, qcell[1], qcell[2], qcell[3])
    ifc2_tensor = Array{Float64}(undef, ifc2_tensor_dims)

    Logging.@debug """
    Allocated the force-constant tensor with following dimensions
    ifc2_tensor_dims = (3, 3, nat, nat, qcell[1], qcell[2], qcell[3])\n
    """ ifc2_tensor_dims

    # Reading the force constant-tensor.
    for i = 1:nfc2
        # Reading the atom displacements and atom/species addresses
        i_displ, j_displ, i_at_addr, j_at_addr =
            parse.(Int64, split(readline(ifc2_file)))

        # Reading the supercell coordinates and actual ifc2-tensor elements
        for j = 1:no_qcells
            temp = split(readline(ifc2_file))

            i_qcell, j_qcell, k_qcell = parse.(Int64, temp[1:(end - 1)])
            ifc2_tensor[
                i_displ,
                j_displ,
                i_at_addr,
                j_at_addr,
                i_qcell,
                j_qcell,
                k_qcell,
            ] = parse(Float64, temp[end])
        end
    end

    close(ifc2_file)

    return ifc2_tensor
end


# End of module Phonon
end
