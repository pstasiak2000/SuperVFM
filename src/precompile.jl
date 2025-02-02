########## Precompilation of some functions for speed

### Straight line vortex using the Schwarz Model
PrecompileTools.@compile_workload begin
    IntPrec = Int32
    FloatPrec = Float32

    α = GetSchwarzTempCoeffs(1.9)

    DimParams = DimensionalParams()
    PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams; FilamentModel=SchwarzModel(α[1],α[2]), IO=IOBuffer())

    Run(PARAMS)
end


### Pre-compile initial conditions not previously tested
PrecompileTools.@compile_workload begin
    IntPrec = Int32
    FloatPrec = Float32 

    initial_conditions = [SingleRing(1.0), SingleHelix(0.1, 0.2, 2π)]

    ### Precompile all the initial vortex conditions
    itC = 0;
    for initf ∈ initial_conditions
        DimParams = DimensionalParams()
        PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams; initf=initf, IO=IOBuffer())

        show(IOBuffer(), PARAMS)

        Run(PARAMS)
    end
end