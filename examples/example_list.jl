
_NAMES_TO_ACCEPTEDNAMES = Dict{Symbol,Array{Symbol}}(
    :FavettoSamson => [:favettosamson, :favetto_samson, :favettoandsamson, :favetto_and_samson],
    :FavettoSamsonAux => [:favettosamsonaux, :favetto_samson_aux, :favettoandsamsonaux, :favetto_and_samson_aux],
    :FitzHughNagumo => [:fitzhughnagumo, :fitz_hugh_nagumo],
    :FitzHughNagumoAlt => [:fitzhughnagumoalt, :fitz_hugh_nagumo_alt],
    :FitzHughNagumoConjug => [:fitzhughnagumoconjug, :fitz_hugh_nagumo_conjug],
    :FitzHughNagumoAux => [:fitzhughnagumoaux, :fitz_hugh_nagumo_aux],
    :FitzHughNagumoAltAuxSimple => [:fitzhughnagumoaltauxsimple, :fitz_hugh_nagumo_alt_aux_simple],
    :FitzHughNagumoAltAuxLin => [Symbol(s, e) for s in [:fitzhughnagumoaltaux, :fitz_hugh_nagumo_alt_aux_] for e in [:lin, :linear, :linearised]],
    :FitzHughNagumoConjugAuxSimple => [:fitzhughnagumoconjugauxsimple, :fitz_hugh_nagumo_conjug_aux_simple],
    :FitzHughNagumoConjugAuxLin => [Symbol(s, e) for s in [:fitzhughnagumoconjugaux, :fitz_hugh_nagumo_conjug_aux_] for e in [:lin, :linear, :linearised]],
    :JansenRit => [:jansenrit, :jansen_rit, :jansen_and_rit, :jansenandrit],
    :JansenRitAux => [:jansenritaux, :jansen_rit_aux, :jansen_and_rit_aux, :jansenandritaux],
    :Lorenz => [:lorenz],
    :LorenzAux => [:lorenzaux, :lorenz_aux],
    :Lorenz96 => [:lorenz96, :lorenz_96],
    :Lorenz96Aux => [:lorenz96aux, :lorenz_96_aux],
    :LotkaVolterra => [:lotka_volterra, :lotkavolterra],
    :LotkaVolterraAux => [:lotka_volterra_aux, :lotkavolterraaux],
    :Prokaryote => [:prokaryote, :prokaryotic],
    :ProkaryoteAux => [:prokaryoteaux, :prokaryoticaux, :prokaryote_aux, :prokaryotic_aux],
    :Sine => [:sine],
    :SineAux => [:sineaux, :sine_aux],
    :SIR => [:sir],
    :SIRAux => [:siraux, :sir_aux],
)

_ADMISSIBLENAMES = collect(Iterators.flatten(values(_NAMES_TO_ACCEPTEDNAMES)))

_ACCEPTEDNAMES_TO_NAMES = Dict{Symbol,Symbol}()
for key in keys(_NAMES_TO_ACCEPTEDNAMES)
    for v in _NAMES_TO_ACCEPTEDNAMES[key]
        _ACCEPTEDNAMES_TO_NAMES[v] = key
    end
end

_NAMES_TO_PATH = Dict{Symbol, String}(
    :FavettoSamson => joinpath("favetto_samson","favetto_samson.jl"),
    :FavettoSamsonAux => joinpath("favetto_samson","favetto_samson_aux.jl"),
    :FitzHughNagumo => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo.jl"),
    :FitzHughNagumoAlt => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt.jl"),
    :FitzHughNagumoConjug => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjugate.jl"),
    :FitzHughNagumoAux => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_aux.jl"),
    :FitzHughNagumoAltAuxSimple => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt_aux_simple.jl"),
    :FitzHughNagumoAltAuxLin => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_alt_aux_linearised.jl"),
    :FitzHughNagumoConjugAuxSimple => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjug_aux_simple.jl"),
    :FitzHughNagumoConjugAuxLin => joinpath("fitz_hugh_nagumo", "fitz_hugh_nagumo_conjug_aux_linearised.jl"),
    :JansenRit => joinpath("jansen_and_rit","jansen_and_rit.jl"),
    :JansenRitAux => joinpath("jansen_and_rit","jansen_and_rit_aux.jl"),
    :Lorenz => joinpath("lorenz_system", "lorenz_system.jl"),
    :LorenzAux => joinpath("lorenz_system","lorenz_system_aux.jl"),
    :Lorenz96 => joinpath("lorenz_96","lorenz_96.jl"),
    :Lorenz96Aux => joinpath("lorenz_96","lorenz_96_aux.jl"),
    :LotkaVolterra => joinpath("lotka_volterra", "lotka_volterra.jl"),
    :LotkaVolterraAux => joinpath("lotka_volterra", "lotka_volterra_aux.jl"),
    :Prokaryote => joinpath("prokaryote","prokaryote.jl"),
    :ProkaryoteAux => joinpath("prokaryote","prokaryote_aux.jl"),
    :Sine => joinpath("sine", "sine.jl"),
    :SineAux => joinpath("sine", "sine_aux.jl"),
    :SIR => joinpath("sir", "sir.jl"),
    :SIRAux => joinpath("sir", "sir_aux.jl"),
)
