# [Sampling of diffusion bridges](@id tutorial_sampling_diff_bridges)
****
> In this tutorial we will illustrate how to sample diffusion bridges using rejection sampling, as well as importance sampling with Brownian bridge proposals. We will illustrate the techniques on the examples of a bivariate double-well potential and the Heston model.

!!! warning
    Sampling of diffusion bridges using rejection sampling or importance sampling with Brownian bridge proposals **is—in general—not an efficient method, and often, it is not even applicable for plenitude of problems**. It is, however, an excellent starting point for learning about simulation of conditioned diffusions. The tutorial is meant as a didactic explanation of the common mechanisms behind most samplers of conditioned diffusions. Nevertheless, it is NOT meant to be used for production purposes. For that we refer the user to [GuidedProposals.jl](https://juliadiffusionbayes.github.io/GuidedProposals.jl/dev/), which has been designed to specifically address the problem of sampling conditioned diffusions.

## Why is there so much text in this tutorial?
---
Sampling of diffusion bridges might seem like an innocent problem if you haven't seen it before. It is however a considerable obstacle that many mathematicians, statisticians and physicists have been attacking for tens of years. And for good reasons—efficient methods that solve this problem are considered to be extremely valuable tools for a variety of real-world applications. No universal solution has been found, however, there are a few contending methodologies that perform particularly well across a variety of problems of this kind. Many methodologies are far too complicated to explain them in a single tutorial. Fortunately, simple versions of rejection and importance sampling fit the bill. Even these however require from us to cover some additional details. Hence the length of the tutorial.

## Rejection sampling
----
Rejection sampling is a simple, but extremely powerful technique for obtaining independent samples from some probability law $\mathbb{P}$ (with density $\dd \mathbb{P}$) using independent samples from another probability law $\mathbb{Q}$ (with density $\dd \mathbb{Q}$). You may, for instance, check out the [wikipedia article](https://en.wikipedia.org/wiki/Rejection_sampling) to learn more.

The basic algorithm is as follows:
```julia
struct RejectionSampling end
# these need to be appropriately defined
function density(P, at=x) ... end
function some_max_bound_on_dP_over_dQ(P, Q) ... end
P, Q = ...

function rand(P, Q, ::RejectionSampling)
    M = some_max_bound_on_dP_over_dQ(P, Q)
    while true
        X° = rand(Q) # proposal
        acceptance_prob = density(P, at=X) / (density(Q, at=X) * M)
        # accept/reject; if reject then re-sample and retry; if accept then return X°
        if rand() ≦ acceptance_prob
            return X°
        end
    end
end
```
A sample returned by the algorithm above:
```julia
X = rand(P, Q, RejectionSampling())
```

is distributed exactly according to $\mathbb{P}$ even though all sampling was done exclusively through the law $\mathbb{Q}$.

!!! tip "Definition"
    The law $\mathbb{P}$ is called the **target** law, whereas $\mathbb{Q}$ is called the **proposal** law.

## Rejection sampling on a path space
---
Rejection sampling on a path space is a fancy name for the algorithm above as applied to diffusion processes, i.e. when the laws $\mathbb{P}$ and $\mathbb{Q}$ are some diffusion laws.

Computation of each individual term in the **acceptance probability**
```julia
acceptance_prob = density(P, at=X) / (density(Q, at=X) * M)
```
which in a mathematical notation becomes:
```math
\frac{1}{M}\frac{\dd \mathbb{P}}{\dd \mathbb{Q}}(X)
```
is not possible when working with diffusion processes. Nevertheless, the overall term can sometimes be found.


!!! note "Some mathematical details"
    In particular, let's assume that the target diffusion law is of the following form:
    ```math
    \begin{equation}\label{eq:target_rejection}
    \dd X_t = α(X_t)\dd t + \dd W_t,\quad X_0=x_0, \quad t\in[0,T],
    \end{equation}
    ```
    with $α:\RR^d→\RR^d$ (at least once continuously differentiable), a volatility coefficient given by the identity matrix and $W$ denoting a $d$-dimensional Brownian motion. Additionally, let's assume that there exists a function $A:\RR^d→\RR$ such that:
    ```math
    ∇A(x) = α(x).
    ```
    (We refer to $A$ as a **potential** function). Then, if we let $\mathbb{P}$ denote the law of \eqref{eq:target_rejection} conditioned on an end-point, i.e. of:
    ```math
    \dd X_t = α(X_t)\dd t + \dd W_t,\quad X_0=x_0,\quad X_T=x_T, \quad t\in[0,T],
    ```
    and we let $\mathbb{Q}$ denote the law of a $d$–dimensional Brownian bridge joining $x_0$ and $x_T$ on $[0,T]$, then the **Radon–Nikodym derivative** between the two is proportional to:
    ```math
    \begin{equation}\label{eq:rnd}
    \frac{\dd \mathbb{P}}{\dd \mathbb{Q}}(X)\propto \exp\left\{- \int_0^Tϕ(X_t)\dd t \right\}\leq 1,
    \end{equation}
    ```
    where
    ```math
    \phi(x):=\frac{1}{2}\left[ \alpha'\alpha + Δ A \right](x) - l.
    ```
    and $l$ is any constant, which we may take in particular to be
    ```math
    l\leq \inf_{x\in\RR^d}\left\{ \frac{1}{2}\left[ \alpha'\alpha + Δ A \right](x) \right\},
    ```
    with the infimum on the right is **assumed to exist**.

The main take-away message from the box above is that if the stochastic differential equation admits a nice enough expression \eqref{eq:target_rejection}, then it is possible to find a closed form formula for the computation of the acceptance probability
```math
\frac{1}{M}\frac{\dd \mathbb{P}}{\dd \mathbb{Q}}(X),
```
and rejection sampling on a path space can be performed by simply repeatedly drawing from $\mathbb{Q}$–the law of Brownian bridges ([see this how-to-guide](@ref how_to_guides) for more info)—and accepting the samples with probability proportional to \eqref{eq:rnd}.

!!! note "Lamperti transformation"
    The restriction to identity matrix volatility coefficient may be slightly relaxed thanks to existence of Lamperti transformation, which allows one to transform a diffusion with a non-identity volatility coefficient to the one with identity volatility. This transform is applicable to all scalar diffusions, however, in multiple dimensions it is more restrictive. See [Ait Sahalia 2008] for more details.

### Example
---
Let's look at a simple example of a bivariate double-well potential model solving the following SDE
```math
\dd X_t = -\frac{1}{2}∇A(X_t)\dd t + \dd W_t,\quad X_0 = x_0,\quad t\in[0,T],
```
where
```math
A(x) := ρ_1\left[ (x_2)^2 - μ_1 \right]^2 + ρ_2\left[ x_2 - μ_2 x_1 \right]^2,
```
and $ρ_1,ρ_2,μ_1,μ_2>0$ are parameters.

```julia

```


Some tedious calculations reveal that:
```julia

```

TBC
