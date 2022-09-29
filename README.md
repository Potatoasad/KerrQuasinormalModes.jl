# KerrQuasinormalModes

This is a package that not only computes Quasinormal Modes for Kerr Black holes, but focuses on an easy to use (yet performant) interface for working with the Quasinormal Mode Functions. These functions are basically solutions to the Tuekolsky Equations for different s,l,m,n and a values. We implement the [Cook-Zalutskiy solver](https://arxiv.org/abs/1410.7698) to compute the quasinormal mode frequency and also the radial expansion coefficients, as well as the spectral coefficients in the angular sector. 

# Usage

Start off by defining a `ModeSequence`: this defines a particular mode without reference to a particular $a$ value. Then you can get the QNM wavefunction corresponding to any particular $a$ value. (In this case $a=0.8$)

```julia
Ω = ModeSequence(s=-2,l=2,m=2,n=0) # Pick a mode
ψ = Ω(0.8) # Get the QNM mode function at a = 0.8
```

## Mode Values

__Values for the Quasinormal Mode Quantities__

```julia
ψ.ω # = 0.5860169749084593 - 0.0756295523345646im
ψ.Alm # = 2.5852939720024275 + 0.20529748567633602im
```

## Wavefunction Evaluation

__QNM Wavefunction values in Boyer Lingquist Coordinates__

```julia
ψ(r) # Returns the value of the radial wavefuntion at r
ψ(r,cos(θ)) # Returns the value of the full wavefunction at r=r, cos(θ)=cos(θ), ϕ=0 and t=0
ψ(r,cos(θ),ϕ) # Returns the value of the full wavefunction at r=r, cos(θ)=cos(θ), ϕ=ϕ and t=0
ψ(r,cos(θ),ϕ,t) # Returns the value of the full wavefunction at r=r, cos(θ)=cos(θ), ϕ=ϕ and t=t
```

__Radial & Spin Weighted Harmonic Wavefunctions__

```julia
ψ.R # Returns the Radial Wavefunction
ψ.R(2.2) # Value of the Radial Wavefunction @ r = 2.2
ψ.R.r₊, ψ.R.r₋  # Radial coordinate of Outer and Inner Horizons
ψ.S # Returns the Spin-Weighted Spheroidal Harmonic Wavefunction
ψ.S(0.7) # Value of the Spin-Weighted Spheroidal Harmonic @ cos(θ) = 0.7
ψ.S(0.7, π/2) # Value of the Spin-Weighted Spheroidal Harmonic @ cos(θ) = 0.7 & ϕ = π/2
ψ.S.Cllʼ # Returns the vector of Spherical-Spheroidal mixing coefficients
```

## Performant Linear Symbolics of QNM Wavefunctions

One can also take _Derivatives_ and _Linear Combinations_ of these modes to emulate the action of differential operators on these functions in a simple way that doesn't come with the slow-down of a general symbolic system. It allows:

- __Derivatives__ of any mode with respect to any of the Boyer Lingquist coordinates `∂r,∂θ,∂ϕ,∂t `
  - Use `∂r,∂θ,∂ϕ,∂t` to create a new Heun-like function. This new function is also very performant during evaluation.
-  __Linear Combinations__ of any mode or it's derivative: (`∂r(Ψ) + m*Ψ  -∂θ(∂θ(Ψ)) `)
  - Doing something like `ψnew = ∂r(Ψ) + m*Ψ  -∂θ(∂θ(Ψ))` is a completely valid operation, and evaluating it as if it was a mode function, like `ψnew(0.6,0.4)`, will be performant.

An example:

```julia
julia> ψnew = ∂r(Ψ) + 3*Ψ  -∂θ(∂θ(Ψ)) # Completely valid expression that is saved as a Linear Combination of Heun-Like functions as shown below
(-1.0 - 0.0im)₋₂Ψ₂₂₀ +
(-1.224744871391589 - 0.0im)₀Ψ₂₂₀ +
(0.2293785997558912 + 1.7983211937744943im)₋₂Ψ₂₂₀ +
(-1.4014125495728098 + 1.352937910894635im)₋₂Ψ₂₂₀ +
(4.075629552334565 + 0.5860169749084593im)₋₂Ψ₂₂₀

julia> ψnew(2.5,0.7)  # Evaluating this new variable is a completely valid and performant operation
0.5945305105574497 + 0.6407097912091125im
```

## Animated Plot of Radial Modes

![RadialModePlot](.docs/QnmAnimated.gif)
