include("./../../src/parsers.jl")
include("./../../src/phonon.jl")

qedata = Parsers.dfpt_qeOutputData("../../examples/espresso.ifc2")
ebdata = Parsers.ebInputData("../../examples/input.nml")

ifc2 = qedata.properties["ifc2"]
