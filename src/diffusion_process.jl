#===============================================================================

            The main macro of this package `@diffusion_process`
            used for defining diffusion processes

===============================================================================#
using MacroTools, StaticArrays

# define some names for the document
_DIMENSIONS = [:dimensions, :dim, :dims, :dimension]
_PARAMETERS = [:parameters, :param, :params]
_CONJUGATE = [:conjugate]
_EXTRA = [:additional, :extra]
_WIENER = [:wiener, :noise, :gaus, :gaussian]
_PROCESS = [:proc, :process, :state, :statespace]
_STATESPACE = [:statespace, :state, :domain, :domains]
_CONSTDIFF = [
    :constdiff,
    :constvola,
    :constdiffusivity,
    :constvolatility,
    :constantdiff,
    :constantvola,
]
_LINEAR = [:linear, :lineardiffusion]
_ELTYPE = [:eltype]
_NUMNONHYPO = [:num_non_hypo, :numnonhypo]
_PHI = [:phi, :ϕ, :φ]

"""
    lowercase(s::Symbol)

Lowercase all letters in a symbol
"""
lowercase(s::Symbol) = Symbol(lowercase(string(s)))

"""
    _symbol_in(s::Symbol, symbols)

Check if symbol `s` is listed in a list of symbols `symbols`
"""
_symbol_in(s::Symbol, symbols) = lowercase(s) in symbols

"""
    _symbol_in(s::QuoteNode, symbols)

Check if quote of a symbol `s` is listed in a list of symbols `symbols`
"""
_symbol_in(s::QuoteNode, symbols) = _symbol_in(eval(s), symbols)

"""
    _symbol_in(::Any, ::Any)

Return `false` by default i.e. if the first argument is not a symbol nor its
quote
"""
_symbol_in(::Any, ::Any) = false

"""
    diffusion_process(name, ex::Expr, p...)

Defines a diffusion process according to a template described in the
documentation of the github repository:
https://github.com/mmider/DiffusionDefinition.jl
"""
macro diffusion_process(name, ex::Expr, p...)
    parse_process(name, MacroTools.striplines(ex), p)
end

"""
    load_diffusion()

Displays available choices of predefined diffusion processes that can be loaded
"""
macro load_diffusion()
    println("The following diffusions have been defined inside the package ",
            "DiffusionDefinition.jl and can be loaded by calling a macro ",
            "@load_diffusion name-of-a-diffusion:")
    for name in keys(_NAMES_TO_ACCEPTEDNAMES)
        println("- $name")
    end
    nothing
end

"""
    load_diffusion(name)

Loads the predefined diffusion process.
"""
macro load_diffusion(name)
    if _symbol_in(name, _ADMISSIBLENAMES)
        if typeof(name) <: QuoteNode
            name = eval(name)
        end
        homedir = joinpath(@__DIR__, "..", "examples")
        importname = _ACCEPTEDNAMES_TO_NAMES[name]
        filepath = _NAMES_TO_PATH[importname]
        path = joinpath(homedir, filepath)
        isfile(path) && include(path)
        !isfile(path) && println("Error, diffusion $name is supposed to be ",
                                "defined but it seems the file does not exist")
        return Meta.parse("import DiffusionDefinition.$importname")
    else
        println("Diffusion $name does not seem to be defined...")
    end
    nothing
end

"""
    parse_process(name , ex::Expr, ::Any)

Parse a template defining a diffusion process, create a corresponding struct and
specified functions, evaluate them in the environment of a package and then
import the struct name to `Main` scope, in which the package has been imported
to.
"""
function parse_process(name , ex::Expr, ::Any)
    p = (
        name = name,
        parameters = Vector{Tuple{Symbol,Union{Symbol,Expr},Symbol,Int64}}(
            undef,
            0
        ),
        dims = Dict{Symbol,Int64}(),
        extras = Dict{Symbol, Any}(),
        fns = Vector{Expr}(undef, 0),
    )
    # parse lines defining parameters to get all parameter names first
    parse_lines!(ex, p, x->(x == _PARAMETERS[1]))
    # parse all other lines
    parse_lines!(ex, p, x->(x != _PARAMETERS[1]))
    fill_unspecified_with_defaults(p)

    abstract_type = prepare_abstract_type(
        p.extras[_LINEAR[1]] ? :LinearDiffusion : :DiffusionProcess,
        p.dims,
        p.extras[_ELTYPE[1]],
        p.extras[_STATESPACE[1]],
    )
    struct_def = createstruct(abstract_type, name, p.parameters)
    add_constdiff_function!(p.fns, p)
    eval(struct_def)
    for fn in p.fns
        eval(fn)
    end
    Meta.parse("import DiffusionDefinition.$name")
end

"""
    parse_lines!(ex::Expr, p, condition)

Parse all lines of the expression `ex`, but process only those which satisfy
`condition`. `p` is a passed-around structure that accumulates processed
information.
"""
function parse_lines!(ex::Expr, p, condition)
    current_label = nothing
    for line in ex.args
        current_label = update_label(line, current_label)
        typeof(line) <: Union{Symbol, QuoteNode} && (continue)
        @assert line.head == :-->
        condition(current_label) && parse_line!(Val{current_label}(), line, p)
    end
end

"""
    update_label(line, current_label)

Update the label, which signifies what type of information a given line in a
template is supposed to be encoding.
"""
function update_label(line, current_label)
    _symbol_in(line, _DIMENSIONS) && return _DIMENSIONS[1]
    _symbol_in(line, _PARAMETERS) && return _PARAMETERS[1]
    _symbol_in(line, _CONJUGATE) && return _CONJUGATE[1]
    _symbol_in(line, _EXTRA) && return _EXTRA[1]
    current_label
end

#------------------------------------------------------------------------------#
#
#                          For parsing: PARAMETERS
#
#------------------------------------------------------------------------------#
"""
    parse_line!(::Val{:parameters}, line, p)

Parse a line that defines parameters of the diffusion. The line must be in a
format:
    name --> parameter-description
"""
function parse_line!(::Val{:parameters}, line, p)
    data_type = line.args[1]

    if typeof(data_type) <: Symbol
        parse_param_single_name(line, p)
    else
        @assert data_type.head == :tuple
        parse_param_multi_names(line, p)
    end
end

"""
    parse_param_single_name(line, p)

Parse a line that defines parameters of the diffusion. The line must be in one
of the formats:
    parameter_name --> (number_of_parameters, datatype)
    parameter_name --> datatype
In the former case defines `number_of_parameters`-many parameters, with names
`parameter_name`i and of `datatype` type. In the latter case defines a
single parameters with name `parameter_name` and of type `datatype`.
"""
function parse_param_single_name(line, p)
    name_stem, disambig_idx, generic_name = get_name_stem(line.args[1], p.parameters)

    data_type = line.args[2]

    if typeof(eval(data_type)) <: DataType
        num_params = 1
    else
        @assert data_type.head == :tuple
        num_params, data_type = data_type.args[1], data_type.args[2]
        generic_name = true
    end

    for i in 1:num_params
        name = ( generic_name ? Symbol(name_stem, disambig_idx+i) : name_stem )
        append!(p.parameters, [(name, data_type, name_stem, disambig_idx+i)])
    end
end

"""
    get_name_stem(name_stem::Symbol, parameters)

Get the stem of a name for a parameters and then add a disambiguation index.
Underscore `_` used in place of name is defaulted to `p`.
"""
function get_name_stem(name_stem::Symbol, parameters)
    if name_stem == :_
        name_stem = :p
    end
    disambig_idx = highest_idx_used(name_stem, parameters)
    name_stem, disambig_idx, (name_stem == :p)
end

"""
    highest_idx_used(name_stem, params)

Find the highest disambiguation index that has been used thus far for a given
`name_stem`.
"""
function highest_idx_used(name_stem, params)
    idxs = [p[4] for p in params if name_stem == p[3]]
    length(idxs) == 0 && return 0
    maximum(idxs)
end

"""
    parse_param_multi_names(line, p)

Parse a line that defines parameters of the diffusion. The line must be in one
of the formats:
    (p_name1, p_name2, ...) --> (number_of_parameters, datatype)
    (p_name1, p_name2, ...) --> datatype
    (p_name1, p_name2, ...) --> (datatype1, datatype2, ...)
In the former two cases defines `number_of_parameters`-many parameters, with
names `p_name1`, `p_name2`, ... and of `datatype` type. In the last case the
datatypes differ from parameter to parameter.
"""
function parse_param_multi_names(line, p)
    names = line.args[1].args
    data_types = line.args[2]

    if typeof(eval(data_types)) <: DataType
        data_types = repeat([data_types], length(names))
    else
        @assert data_types.head == :tuple
        if typeof(eval(data_types.args[1])) <: Number
            @assert length(data_types.args) == 2
            @assert eval(data_types.args[1]) == length(names)
            @assert typeof(eval(data_types.args[2])) <: DataType
            data_types = repeat([data_types.args[2]], length(names))
        else
            @assert length(data_types.args) == length(names)
            @assert all(map(x->(typeof(eval(x))<:DataType), data_types.args))
            data_types = data_types.args
        end
    end

    current_param_names = [param[1] for param in p.parameters]
    for (n,d) in zip(names, data_types)
        @assert !(n in current_param_names)
        append!(p.parameters, [(n,d,n,1)])
    end
end

#------------------------------------------------------------------------------#
#
#                       For parsing: CONJUGATE UPDATES
#
#------------------------------------------------------------------------------#

"""
    parse_line!(::Val{:conjugate}, line, p)

Parse a line that defines conjugate updates.
"""
function parse_line!(::Val{:conjugate}, line, p)
    if typeof(line.args[1]) <: Expr
        @assert line.args[1].head == :call
        function_name = line.args[1].args[1]
        fns = Vector{Expr}(undef, 0)
        if _symbol_in(function_name, _PHI)
            t_sym = line.args[1].args[2]
            x_sym = line.args[1].args[3]
            @assert line.args[2].head == :tuple
            num_funs = length(line.args[2].args)
            for expr in line.args[2].args
                add_phi_function!(expr, fns, t_sym, x_sym, p)
            end
        elseif function_name == :nonhypo
            add_nonhypo_function!(fns, line, p)
        end
        append!(p.fns, fns)
    elseif _symbol_in(line.args[1], _NUMNONHYPO)
        num = line.args[2]
        @assert typeof(num) <: Integer
        p.extras[_NUMNONHYPO[1]] = num
    end
end

"""
    add_nonhypo_function!(fns, line, p)

Add a definition of a function `nonhypo` that takes a current diffusion state
`x` and returns the elliptic coordinates.
"""
function add_nonhypo_function!(fns, line, p)
    fn_def = Expr(:call,
        :nonhypo,
        Expr(:(::),
            :P,
            p.name
        ),
        line.args[1].args[2]
    )
    push!(fns, Expr(:(=), fn_def, line.args[2]))
end

"""
    add_phi_function!(expr::Expr, fns::Vector{Expr}, t::Symbol, x::Symbol, p)

Construct a `phi` function used for conjugate updates.
"""
function add_phi_function!(
    expr::Expr,
    fns::Vector{Expr},
    t::Symbol,
    x::Symbol,
    p
    )
    idx = length(fns)
    name = p.name
    phi_def = Expr(:call,
        _PHI[1],
        Expr(:(::),
            Val{idx}
        ),
        :t,
        :x,
        Expr(:(::), :P, p.name)
    )

    @assert expr.head == :tuple
    @assert length(expr.args) == p.dims[_PROCESS[1]]
    phi_body = copy(expr)
    for i in 1:p.dims[_PROCESS[1]]
        cleanup_param_names!(phi_body.args[i], t, x, map(x->x[1], p.parameters))
    end

    phi = Expr(:(=), phi_def, phi_body)
    append!(fns, [phi])
end

"""
    cleanup_param_names!(expr, t::Symbol, x::Symbol, params)

A hepler function for cleaning up the names of variables in the definition of a
`phi` function.
"""
function cleanup_param_names!(expr, t::Symbol, x::Symbol, params)
    for (i,arg) in enumerate(expr.args)
        if typeof(arg) <: Symbol
            arg == t && (expr.args[i] = :t)
            arg == x && (expr.args[i] = :x)
            arg in params && (expr.args[i] = Expr(:., :P, QuoteNode(arg)))
        elseif typeof(arg) == Expr
            cleanup_param_names!(arg, t, x, params)
        end
    end
end

#------------------------------------------------------------------------------#
#
#                       For parsing: ADDITIONAL INFO
#
#------------------------------------------------------------------------------#

"""
    parse_line!(::Val{:additional}, line, p)

Parse a line that defines additional information about a diffusion process.
"""
function parse_line!(::Val{:additional}, line, p)
    name = line.args[1]
    _symbol_in(name, _STATESPACE) && (name = _STATESPACE[1])
    _symbol_in(name, _CONSTDIFF) && (name = _CONSTDIFF[1])
    _symbol_in(name, _LINEAR) && (name = _LINEAR[1])
    _symbol_in(name, _ELTYPE) && (name = _ELTYPE[1])

    p.extras[name] = line.args[2]
end

"""
    add_constdiff_function!(fns, p)

Add a definition of a function `consdiff` that indicates if the diffusion
coefficient is constant
"""
function add_constdiff_function!(fns, p)
    fn_def = Expr(:call,
        :constdiff,
        Expr(:(::),
            :P,
            p.name
        )
    )
    push!(fns, Expr(:(=), fn_def, p.extras[_CONSTDIFF[1]]))
end

#------------------------------------------------------------------------------#
#
#                   For parsing: DIFFUSION DIMENSION
#
#------------------------------------------------------------------------------#

"""
    parse_line!(::Val{:dimensions}, line, p)

Parse a line that defines the dimension of a diffusion process and the driving
Brownian motion.
"""
function parse_line!(::Val{:dimensions}, line, p)
    dim_of_what, dim = line.args[1], line.args[2]
    _symbol_in(dim_of_what, _WIENER) && (p.dims[_WIENER[1]] = dim)
    _symbol_in(dim_of_what, _PROCESS) && (p.dims[_PROCESS[1]] = dim)
end

#------------------------------------------------------------------------------#
#
#                            final routines
#
#------------------------------------------------------------------------------#

"""
    fill_unspecified_with_defaults(p)

Fill all unspecified variables with default values
"""
function fill_unspecified_with_defaults(p)
    fill_unspecified_with_defaults(Val{:dimensions}(), p)
    fill_unspecified_with_defaults(Val{:additional}(), p)
end

"""
    fill_unspecified_with_defaults(::Val{:dimensions}, p)

If unspecified, the dimension of the stochastic process and the driving Brownian
motion is set to 1.
"""
function fill_unspecified_with_defaults(::Val{:dimensions}, p)
    for key in [_WIENER[1], _PROCESS[1]]
        !haskey(p.dims, key) && (p.dims[key] = 1)
    end
end

"""
    fill_unspecified_with_defaults(::Val{:additional}, p)

If unspecified, there are no restriction on a state space, the volatility
coefficient is assumed constant, the diffusion is not linear and the datatype
of each coordinate is set to `Float64`.
"""
function fill_unspecified_with_defaults(::Val{:additional}, p)
    defaults = [
        (_STATESPACE[1], UnboundedStateSpace()),
        (_CONSTDIFF[1], false),
        (_LINEAR[1], false),
        (_ELTYPE[1], Float64),
        (_NUMNONHYPO[1], (
            haskey(p.extras, _STATESPACE[1]) ?
            p.extras[_STATESPACE[1]] :
            1
        )),
    ]
    for (key, default_val) in defaults
        !haskey(p.extras, key) && (p.extras[key] = default_val)
    end
end

"""
    prepare_abstract_type(stem, dims, data_type, state_restr)

Create a string defining a parent, abstract type from its `stem`, the dimensions
`dims` of the process and the driving Brownian motion, the datatype `data_type`
of each coordinate and the restrictions on the state space `state_restr`.
"""
function prepare_abstract_type(stem, dims, data_type, state_restr)
    wiener_dim, proc_dim = dims[_WIENER[1]], dims[_PROCESS[1]]
    "$stem{$data_type,$proc_dim,$wiener_dim,$state_restr}"
end

"""
    createstruct(abstract_type, name, params)

Create code that defines a struct defining a diffusion process.
"""
function createstruct(abstract_type, name, params)
    header = "struct $name <: $abstract_type\n"
    params_vec = map(params) do (name, data_type, _, _)
        "\t$name::$data_type\n"
    end
    Meta.parse(string(header, join(params_vec), "end\n"))
end
