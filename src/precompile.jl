PrecompileTools.@compile_workload begin
    IntPrec = Int32
    FloatPrec = Float32

    initial_conditions = [SingleLine(),SingleRing(1.0),SingleHelix(0.1,0.2,2π)]

    for initf ∈ initial_conditions
        DimParams = DimensionalParams()
        PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams)

        Run(PARAMS)
    end
end