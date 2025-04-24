# Notes on the `kickstart.jl` script



## Diagonalization

### CPU

- Using the base implementation from `LinearAlgebra`


### GPU

- `CUDA.CUSOLVER.heevjBatched!` is **H**ermitian Matrix **J**acobi based **E**igen **V**alue **D**ecomposition which comes from the `CUSOLVER` library [here](https://docs.nvidia.com/cuda/cusolver/index.html#cusolverdn-t-syevjbatched)

    > [!CAUTION]
    > What is batching?

## Questions

**1.** Is this bad coding practice because it does not tell us that it mutates
    it's arguments?

~~~julia
function diag_dynmats(dynmats) 
    eigen(dynmats)
    return 1
end
~~~

**2.** What does and expression of the type `CUDA.@sync` mean in front of a variable?

**3.** What is batching over an index?
