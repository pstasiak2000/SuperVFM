var documenterSearchIndex = {"docs":
[{"location":"misc/#Miscelleneus-functions","page":"Miscelleneus functions","title":"Miscelleneus functions","text":"","category":"section"},{"location":"misc/","page":"Miscelleneus functions","title":"Miscelleneus functions","text":"<!– @docs convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments) –>","category":"page"},{"location":"vortex_configs/single_ring/#Single-vortex-ring","page":"Single vortex ring","title":"Single vortex ring","text":"","category":"section"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"Initialises a single, vortex ring of size R.","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"(Image:  Single vortex ring)","category":"page"},{"location":"vortex_configs/single_ring/#Usage","page":"Single vortex ring","title":"Usage","text":"","category":"section"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"struct SingleRing{A} <: InitCond\n    Radius::A\nend","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"Vortex ring in the yz-plane, propagating in the x direction. Use the field Radius to set the radius of the vortex ring.","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"#Set an initial condition of ring size 1.0\nIC = VortexFilament.SingleRing(1.0f0)","category":"page"},{"location":"manual/#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Here everything will be listed, the routines for computing ","category":"page"},{"location":"manual/#Vortex-filament-structure-on-GPU","page":"Manual","title":"Vortex filament structure on GPU","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"In the standard implementations of the Vortex Filament Method (VFM), filaments are constructed by a Fortran 'type' array. Unfortunately, though structure types are native to Julia, structures defined and adapted using the Adapt package are required to be immutable. The implementation of the VFM requires structures to mutate according to the number of vortex points, the value will fluctuate to maintain discretisation between filaments. ","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"tip: Tip\nThe GPU array can be converted to a filament structure by calling convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments). [WIP!!!] The GPU array can be converted to a filament structure by calling convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments).","category":"page"},{"location":"manual/#Float32-Array","page":"Manual","title":"Float32 Array","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Rows Variable Description\n1-3 xyz position of the vortex point\n4-6 u velocity of the vortex point\n7-9 u1 stored velocity for Adams-Bashforth\n10-12 u2 stored velocity for Adams-Bashforth\n13-15 f_u_mf velocity due to mutual friction\n16-18 f_u_sup velocity due to background flow\n19-21 f_u_loc velocity due to local contribution\n22-24 f_u_pll velocity due to flow across tangent\n25-27 f_f_mf mutual friction of the vortex point\n28-30 f_u_n normal fluid velocity\n31 f_curv curvature of vortex point\n32 f_stretch helicity density mathbfscdotmathbfv_n(mathbfs)\n33-35 f_ghosti ghost point infront\n36-38 f_ghostb ghost point behind\n39 f_closestd distance to closest filament","category":"page"},{"location":"manual/#Int32-Array","page":"Manual","title":"Int32 Array","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Rows Variable Description\n1 f_infront vortex point infront\n2 f_behind vortex point behind\n3 f_closest vortex point closest","category":"page"},{"location":"manual/#Variables","page":"Manual","title":"Variables","text":"","category":"section"},{"location":"manual/#Vortex-Configurations","page":"Manual","title":"Vortex Configurations","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Here is a list of all of the vortex configurations that can be used for the vortex filament method. Specific parameters relating to the initialisation of the vortex line are stored in the corresponding structure blocks.","category":"page"},{"location":"manual/#List-of-implemented-initial-conditions","page":"Manual","title":"List of implemented initial conditions","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Pages = Main.SUBSECTION_PAGES\nDepth=1","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"int x dx","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"a = 1\nb = 2\na + b","category":"page"},{"location":"#VortexFilament.jl","page":"Home","title":"VortexFilament.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A little bit of text here to see the changes.","category":"page"},{"location":"#Outline","page":"Home","title":"Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"vortex_configs/straight_line/#Straight-line","page":"Straight line","title":"Straight line","text":"","category":"section"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"Initialises a single, straight line vortex at the centre of the computational box.","category":"page"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"(Image:  Straight line vortex)","category":"page"},{"location":"vortex_configs/straight_line/#Usage","page":"Straight line","title":"Usage","text":"","category":"section"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"struct SingleLine <: InitCond\nend\nAdapt.@adapt_structure SingleLine","category":"page"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"To set up a single vortex line, simply pass SingleLine().","category":"page"}]
}
