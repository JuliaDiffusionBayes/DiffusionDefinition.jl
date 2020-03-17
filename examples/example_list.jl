
_NAMES_TO_ACCEPTEDNAMES = Dict{Symbol,Array{Symbol}}(
    :Lorenz => [:lorenz],
    :LotkaVolterra => [:lotka_volterra, :lotkavolterra],
    :LotkaVolterraAux => [:lotka_volterra_aux, :lotkavolterraaux],
    :FitzHughNagumo => [:fitzhughnagumo, :fitz_hugh_nagumo],
    :FitzHughNagumoAlt => [:fitzhughnagumoalt, :fitz_hugh_nagumo_alt],
    :FitzHughNagumoConjug => [:fitzhughnagumoconjug, :fitz_hugh_nagumo_conjug],
    :FitzHughNagumoAux => [:fitzhughnagumoaux, :fitz_hugh_nagumo_aux],
    :FitzHughNagumoAltAuxSimple => [:fitzhughnagumoaltauxsimple, :fitz_hugh_nagumo_alt_aux_simple],
    :FitzHughNagumoAltAuxLin => [Symbol(s, e) for s in [:fitzhughnagumoaltaux, :fitz_hugh_nagumo_alt_aux_] for e in [:lin, :linear, :linearised]],
    :FitzHughNagumoConjugAuxSimple => [:fitzhughnagumoconjugauxsimple, :fitz_hugh_nagumo_conjug_aux_simple],
    :FitzHughNagumoConjugAuxLin => [Symbol(s, e) for s in [:fitzhughnagumoconjugaux, :fitz_hugh_nagumo_conjug_aux_] for e in [:lin, :linear, :linearised]],
    :Sine => [:sine],
    :SineAux => [:sineaux, :sine_aux],
    :SIR => [:sir],
    :SIRAux => [:sir, :sir_aux],
)

_ADMISSIBLENAMES = collect(Iterators.flatten(values(_NAMES_TO_ACCEPTEDNAMES)))

_ACCEPTEDNAMES_TO_NAMES = Dict{Symbol,Symbol}()
for key in keys(_NAMES_TO_ACCEPTEDNAMES)
    for v in _NAMES_TO_ACCEPTEDNAMES[key]
        _ACCEPTEDNAMES_TO_NAMES[v] = key
    end
end

_NAMES_TO_PATH = Dict{Symbol, String}(
    :Lorenz => joinpath("lorenz_system", "lorenz_system.jl"),
    :LotkaVolterra => joinpath("lotka_volterra", "lotka_volterra.jl"),
    :LotkaVolterraAux => joinpath("lotka_volterra", "lotka_volterra_aux.jl"),
    :FitzHughNagumo => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo.jl"),
    :FitzHughNagumoAlt => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt.jl"),
    :FitzHughNagumoConjug => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjug.jl"),
    :FitzHughNagumoAux => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_aux.jl"),
    :FitzHughNagumoAltAuxSimple => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt_aux_simple.jl"),
    :FitzHughNagumoAltAuxLin => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt_aux_linearised.jl"),
    :FitzHughNagumoConjugAuxSimple => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjug_aux_simple.jl"),
    :FitzHughNagumoConjugAuxLin => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjug_aux_linearised.jl"),
    :Sine => joinpath("sine", "sine.jl"),
    :SineAux => joinpath("sine", "sine_aux.jl"),
    :SIR => joinpath("sir", "sir.jl"),
    :SIRAux => joinpath("sir", "sir_aux.jl"),
)
