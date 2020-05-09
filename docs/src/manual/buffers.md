# [Buffers](@id explain_buffers)
Buffers are simply structs gathering containers that are used for in-place computations. They inherit from
```@docs
DiffusionDefinition.AbstractBuffer
```
and don't require any special methods, except for functions `b!` and `σ!` (or for linear diffusions `B!` and `β!` in place of `b!`) having to specify explicitly how the drift and volatility are supposed to be saved inside the buffers. We provide implementation of standard buffers that should be sufficient for many problems:
```@docs
DiffusionDefinition.StandardEulerBuffer
DiffusionDefinition.LinearDiffBuffer
```
!!! tip
    If you call in-place methods multiple times, then make sure to pre-define a buffer and then pass it explicitly to `solve!`.
