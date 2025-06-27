module Phonon

import Logging

include("./parsers.jl")

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

    for i = 1:3
        for j = 1:3
            for iat = 1:nat
                full_sum = sum(ifc2_tensor[i, j, iat, :, :, :, :])
                ifc2_tensor[i, j, iat, iat, 1, 1, 1] -=
                    ifc2_tensor[i, j, iat, iat, 1, 1, 1] - full_sum
            end
        end
    end
    return
end

"""
From this struct we will draw all the data we need to perform a calculation
"""
struct HarmonicScaffolding
    sqrtmasses_mat::Array{Float64}
    # Square matrix prefactor to the dynamical matrix. Diagonal holds mass 
    # of atom i. All other elements hold mass factors for atoms i and j
    weights::Array{Float64}

    function HarmonicScaffolding(
        ebdata::Parsers.ebInputData,
        dfptdata::Parsers.dfpt_qeOutputData,
    ) end
end

# End of module Phonon
end
