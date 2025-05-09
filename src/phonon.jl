module Phonon

import Statistics

"""
    read_ifc2(ifc2_path::String, ebinput_path::String)

Read a Quantum-Espresso-ifc2-output file into a force constant tensor.

The function needs to read a corresponding elphbolt-input file for crystals
information as well.
"""
function read_ifc2(ifc2_path::String, ebinput_path::String)
    # TODO Check that files exist

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

    # TODO Deal with the special ibrav=0 case

    # NEXT ntype LINE(S)
    mass = Vector{Float64}(undef, ntype)
    label = Array{String}(undef, ntype)
    for i = 1:ntype
        temp = split(readline(ifc2_file))
        println(temp)
        label[i] = chop(temp[2], head = 1, tail = 1)
        mass[i] = parse(Float64, temp[3])
        println(label)
        println(mass)
    end

    # NEXT nat LINE(S)
    # From these lines we get which species (tipo) sits on which basis 
    # position (basis_pos)
    tipo = Vector{Int16}(undef, nat)
    basis_pos = Matrix{Float64}(undef, (nat, 3))
    for i = 1:nat
        temp = split(readline(ifc2_file))
        println(temp)
        tipo[i] = parse(Int16, temp[2])
        basis_pos[i, :] = parse.(Float64, temp[3:end])
    end
    println(basis_pos)
    println(tipo)

    # NEXT LINE
    # Is the material a polar one?
    polar_key = split(readline(ifc2_file))[1]
    println(polar_key)
    if polar_key == "T"
        # NEXT 3 LINES
        # Read the dielectric tensor
        dielectric_tensor = Matrix{Float64}(undef, (3, 3))
        for i = 1:3
            dielectric_tensor[i, :] =
                parse.(Float64, split(readline(ifc2_file)))
        end
        println(dielectric_tensor)

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
    end

    # NEXT LINE
    # Read the super cell sizes
    qcell = parse.(Int32, split(readline(ifc2_file)))
    println(qcell)

    # NEXT 3*3*nat*nat*( qcell[1]*qcell[2]*qcell[3] + 1 ) lines
    # First we define the number of force-constant-supercell-fields 
    nfc2 = 3 * 3 * nat * nat
    no_qcells = qcell[1] * qcell[2] * qcell[3]

    # Allocating the force tensor
    ifc2_tensor =
        Array{Float64}(undef, (3, 3, nat, nat, qcell[1], qcell[2], qcell[3]))

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
            ] = parse(Float64, temp[4])
        end
    end

    close(ifc2_file)

    # Verbose output
    # for i = 1:size(zeff)[1]
    #     println("zeff[$(i)]=", i)
    #     for j = 1:size(zeff)[2]
    #         println(zeff[i, j, :])
    #     end
    # end
end

end
