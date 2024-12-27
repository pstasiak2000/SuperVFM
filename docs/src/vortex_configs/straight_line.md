# Straight line 

*Initialises a single, straight line vortex at the centre of the computational box.*

![ Straight line vortex](../assets/straight_line.svg)

## Usage
```julia
struct SingleLine <: InitCond
end
<<<<<<< HEAD
=======
Adapt.@adapt_structure SingleLine
>>>>>>> b28cc71f9ab0bd7658202bc3e7a22387ef04881e
```

To set up a single vortex line, simply pass `SingleLine()`.