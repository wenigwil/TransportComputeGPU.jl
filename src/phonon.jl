module Phonon

import Logging

include("./parsers.jl")

"""
From this struct we will draw all the data we need to perform a calculation
"""
struct HarmonicScaffolding
    mass_prefactor::Array{Float64}

    # function HarmonicScaffolding(
    #     ebdata::Parsers.ebInputData,
    #     dfptdata::Parsers.dfpt_qeOutputData,
    # ) end
end

"""
    enforce_acoustic_sum_rule!(ifc2_tensor)

Enforce the acoustic sum rule on an tensor of force constants. It can be derived
from Newtons 3rd Law for the atoms in the unitcell at the origin.

The force constants should have the shape `(3, 3, nat, nat, sc[1], sc[2], sc[3])`,
where `nat` are the number of atoms in a unitcell and `sc` is a vector of integers
that correspond to the number of unitcells that make up the supercell over which
the tensor is defined.

# References

  - G. J. Ackland et. al. 1997 "Practical methods in ab initio lattice dynamics"
"""
function enforce_acoustic_sum_rule!(ifc2_tensor::Array{Float64})
    # Grab the number of atoms from the shape
    nat = size(ifc2_tensor)[3]

    for i in 1:3
        for j in 1:3
            for iat in 1:nat
                full_sum = sum(ifc2_tensor[i, j, iat, :, :, :, :])
                ifc2_tensor[i, j, iat, iat, 1, 1, 1] =
                    ifc2_tensor[i, j, iat, iat, 1, 1, 1] - full_sum
            end
        end
    end
    return
end

"""
    build_mass_prefactor(
            numbasisatoms::Int64,
            basisatom2species::Vector{Int64},
            species2mass::Vector{Float64})

Build the mass prefactor that is needed in the computation of the
dynamical matrix.

The mass_prefactor is a matrix of size `(numbasisatoms,  numbasisatoms)`. Its `i`-th diagonal element is just the mass of
`species[basisatom2species[i]]` which we will call `mass[i]` in this
comment. Every other `(i,j)`-th element is `sqrt( mass[i] * mass[j] )`
"""
function build_mass_prefactor(
    numbasisatoms::Int64,
    basisatom2species::Vector{Int64},
    species2mass::Vector{Float64},
)

    # Build a vector that holds the mass for each basisatom
    mass = Vector{Float64}(undef, (numbasisatoms))
    for basisatom in 1:numbasisatoms
        # Grab the species of each basis atom
        basisspecies = basisatom2species[basisatom]
        # Convert the species of each basis atom to its mass
        mass[basisatom] = species2mass[basisspecies]
    end

    mass_prefactor = sqrt.(mass * transpose(mass))

    return mass_prefactor
end

"""
    build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

Build all vectors that connect each position defined in `basis`.

The tensor that is built is anti-symmetric in the first two indices and
thus contains zero-vectors this diagonal. The `(i,j)`-th element of
`basisconnectors` contains the difference of `basis[i]` and `basis[j]` with
the `(j,i)`-th element being the negative of the former mentioned.
"""
function build_basisconnectors(numbasisatoms::Int64, basis::Matrix{Float64})

    # Build square (numatoms x numatoms)-matrix of vectors such that every 
    # row consists of all basis vectors of all basisatoms
    basis_dublicate = Array{Float64}(undef, (numbasisatoms, numbasisatoms, 3))

    # Fill the first row
    basis_dublicate[1, :, :] = basis[:, :]
    # Dublicate the first row to the other
    for j in 1:numbasisatoms
        basis_dublicate[j, :, :] = basis_dublicate[1, :, :]
    end

    # Connectors are the difference between the dublicate and the 
    # "transpose" (as for a matrix vector-elements) of the dublicate 
    basisconnectors = basis_dublicate - permutedims(basis_dublicate, (2, 1, 3))

    return basisconnectors
end

"""
    build_supercell_positions(
        unit_multiplicity_super::Vector{Int64},
        super_multiplicity_ultra::Vector{Int64},
        unit_lattvecs::Matrix{Float64})

Generate an ultracell by repetion of supercells, for which you build
all positions to as well as the positions half length squared.

The supercell is made up of a multiplicity of unitcells in each direction
of `unit_lattvecs`. Analogously the ultracell is made up of a
multiplicity of supercells in each direction of `unit_lattvecs`.

# Important

  - For now `super_multiplicity_ultra` must an odd number
  - It is assumed that the `unit_lattvecs` column (2nd index) goes through the
    cartesian coordinates!
  - It is assumed that the orderings in the `unit_multiplicity_super` and
    `unit_lattvecs` align
"""
function build_supercell_positions(
    unit_multiplicity_super::Vector{Int64},
    super_multiplicity_ultra::Vector{Int64},
    unit_lattvecs::Matrix{Float64},
)

    # Supercell multiplicity must be odd in every direction
    for i in 1:3
        if iseven(super_multiplicity_ultra[i]) || super_multiplicity_ultra[i] < 3
            Logging.@error "Multiplicity of the supercells in the ultracell is not valid!"
            return
        end
    end

    # Build the lattice vectors of the supercell by stretching the unitcell 
    # lattvecs
    super_lattvecs = Matrix{Float64}(undef, (3, 3))
    for row in 1:3
        super_lattvecs[row, :] = unit_lattvecs[row, :] * unit_multiplicity_super[row]
    end

    # The supercell_positions array will not contain the origin (0,0,0)
    num_supercell_positions = prod(super_multiplicity_ultra) - 1
    # Saving the supercell positions in the ultracell
    super_positions = Matrix{Float64}(undef, (num_supercell_positions, 3))
    # Saving the distance of the bisector between every supercell position 
    # and the origin
    super_bisector_dist = Vector{Float64}(undef, num_supercell_positions)

    # In Quantum Espresso ultra_range_max is fixed to 2
    ultra_range_max = div.(super_multiplicity_ultra, 2)
    ultra_range = range.(-ultra_range_max, ultra_range_max)
    # Building the ultracell as a cube with the same number of 
    # supercell_positions in every direction around the origin.
    j = 1
    for m1 in ultra_range[1]
        for m2 in ultra_range[2]
            for m3 in ultra_range[3]
                # Exclude the origin as mentioned before
                if iszero([m1, m2, m3])
                    continue
                end

                super_positions[j, :] = begin
                    super_lattvecs[1, :] * m1 +
                    super_lattvecs[2, :] * m2 +
                    super_lattvecs[3, :] * m3
                end

                super_bisector_dist[j] = begin
                    0.5 * transpose(super_positions[j, :]) * super_positions[j, :]
                end

                j += 1
            end
        end
    end

    return super_positions, super_bisector_dist
end

# End of module Phonon
end
