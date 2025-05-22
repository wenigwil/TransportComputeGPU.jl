import Logging

include("../../src/parsers.jl")

# eb_namelists = Parsers.ebInputData("../../examples/input.nml", 3, 4)
#
# for (k, v) in eb_namelists.allocations
#     println(k, " = ", v)
# end
# for (k, v) in eb_namelists.crystal_info
#     println(k, " = ", v)
# end

data = Parsers.ebInputData("./../../examples/input.nml")


