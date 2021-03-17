# KerrQuasinormalModes

This is a package that not only computes Quasinormal Modes for Kerr Black holes, but focuses on an easy to use (yet performant) interface for working with the Quasinormal Mode Functions. These functions are basically solutions to the Tuekolsky Equations for different s,l,m,n and a values. We implement the [Cook-Zalutskiy solver](https://arxiv.org/abs/1410.7698) to compute the quasinormal mode frequency and also the radial expansion coefficients, as well as the spectral coefficients in the angular sector. 

## Usage
Define a new quasinormal mode using the `qnmfunction` command. 
```julia
Ψ = qnmfunction(s=-2,l=2,m=1,n=1,a=0.2)
```
This computes the quasinormal mode frequency internally and also has built in formulae for the mode function expansions.
```julia
r = 2.0:0.01:8.0      #Sample points from the r = 2M out till r = 8M
Ψ.(r)                 #Radial mode sampled at many r points
```

## Animated Plot of Radial Modes
![RadialModePlot](docs/QnmAnimated.gif)
