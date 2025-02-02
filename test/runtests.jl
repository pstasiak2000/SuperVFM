push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Unitful
using Printf

make_animation = false
make_plot = false

#Vortex initial condition
initf = SingleRing()
# IC = SingleRing(0.25)

IntPrec = Int32
FloatPrec = Float32

DimParams = DimensionalParams()
PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams; initf=initf, IO=IOBuffer())

Run(PARAMS)