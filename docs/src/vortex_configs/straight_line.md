# Straight line 

*Initialises a single, straight line vortex at the centre of the computational box.*

![ Straight line vortex](../assets/straight_line.svg)

## Usage
```julia
struct SingleLine <: InitCond
end
Adapt.@adapt_structure SingleLine
```

To set up a single vortex line, simply pass `SingleLine()`.