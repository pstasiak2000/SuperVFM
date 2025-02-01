### Precompilation of some functions for speed
PrecompileTools.@compile_workload begin
    IntPrec = Int32
    FloatPrec = Float32

    initial_conditions = [SingleLine(),SingleRing(1.0),SingleHelix(0.1,0.2,2π)]

    ### Precompile all the initial vortex conditions
    for initf ∈ initial_conditions
        DimParams = DimensionalParams()
        PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams; initf=initf, IO=IOBuffer())

        show(IOBuffer(),PARAMS)

        Run(PARAMS)
    end
end