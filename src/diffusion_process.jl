#===============================================================================

            The main macro of this package `@diffusion_process`
            used for defining diffusion processes

===============================================================================#
# define some names for the document
_DIMENSIONS = [:dimensions, :dim, :dims, :dimension]
_PARAMETERS = [:parameters, :param, :params]
_CONSTPARAMETERS = [
    :const_parameters, :const_param, :const_params,
    :constparameters, :constparam, :constparams,
    :constant_parameters, :constant_param, :constant_params,
    :constantparameters, :constantparam, :constantparams
]
_CONJUGATE = [:conjugate]
_ENDPOINTINFO = [:aux_info, :auxiliary_info, :end_points, :end_point_info]
_ENDPOINTEXTRA = [
    [:t0, :t_0],
    [:T],
    [:x0, :y0, :state0, :x_0, :y_0, :state_0],
    [:v0, :obs0, :v_0, :obs_0],
    [:vT, :obsT,  :v_T, :obs_T],
    [:xT, :yT, :stateT, :x_T, :y_T, :state_T],
]
_EXTRA = [:additional, :extra]
_WIENER = [:wiener, :noise, :gaus, :gaussian]
_PROCESS = [:proc, :process, :state, :statespace]
_STATESPACE = [:statespace, :state, :domain, :domains]
_CONSTDIFF = [
    :constdiff,
    :constvola,
    :constdiffusivity,
    :constvolatility,
    :constσ,
    :constantdiff,
    :constantvola,
    :constantdiffusivity,
    :constantvolatility,
    :constantσ,
]
_DIAGONALDIFF = [
    :diagonaldiff,
    :diagonalvola,
    :diagonaldiffusivity,
    :diagonalvolatility,
    :diagonalσ,
    :diagdiff,
    :diagvola,
    :diagdiffusivity,
    :diagvolatility,
    :diagσ,
]
_SPARSEDIFF = [
    :sparsediff,
    :sparsevola,
    :sparsediffusivity,
    :sparsevolatility,
    :sparseσ,
]
_DIAGONALBMAT = [
    :diagonalbmat,
    :diagonalb,
    :diagonalbmatrix,
]
_SPARSEBMAT = [
    :sparsebmat,
    :sparseb,
    :sparsebmatrix,
]
_LINEAR = [:linear, :lineardiffusion]
_ELTYPE = [:eltype]
_NUMNONHYPO = [:num_non_hypo, :numnonhypo]
_PHI = [:phi, :ϕ, :φ]

_TEMP_DIMENSION_PROCESS = -1
_TEMP_DIMENSION_WIENER = -1

"""
    Base.lowercase(s::Symbol)

Lowercase all letters in a symbol
"""
Base.lowercase(s::Symbol) = Symbol(lowercase(string(s)))

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
    found_me, importname = _load_diffusion(name)
    found_me && return Meta.parse("import DiffusionDefinition.$importname")
    nothing
end

function _load_diffusion(name)
    if _symbol_in(name, _ADMISSIBLENAMES)
        if typeof(name) <: QuoteNode
            name = eval(name)
        end
        homedir = joinpath(@__DIR__, "..", "examples")
        importname = _ACCEPTEDNAMES_TO_NAMES[lowercase(name)]
        filepath = _NAMES_TO_PATH[importname]
        path = joinpath(homedir, filepath)
        isfile(path) && include(path)
        !isfile(path) && println("Error, diffusion $name is supposed to be ",
                                "defined but it seems the file does not exist")
        return true, importname
    end
    println("Diffusion $name does not seem to be defined...")
    false, nothing
end

macro load_variable_diffusion(name, dim_proc, args...)
    global _TEMP_DIMENSION_PROCESS = eval(dim_proc)
    global _TEMP_DIMENSION_WIENER = (
        length(args) > 0 ?
        eval(args[1]) :
        _TEMP_DIMENSION_PROCESS
    )
    found_me, importname = _load_diffusion(name)
    found_me && return Meta.parse("import DiffusionDefinition.$importname")
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
    name, curly_name, template_args = process_name(name)
    p = (
        name = name,
        curly_name = curly_name,
        template_args = template_args,
        parameters = Vector{Tuple{Symbol,Union{Symbol,Expr},Symbol,Int64}}(
            undef,
            0
        ),
        dims = Dict{Symbol,Int64}(),
        extras = Dict{Symbol, Any}(),
        fns = Vector{Expr}(undef, 0),
        constp = Symbol[],
        nonconstp = Symbol[],
    )

    # parse lines defining parameters to get all parameter names first
    parse_lines!(ex, p, x->( x in [_PARAMETERS[1], _CONSTPARAMETERS[1]] ))
    # parse all other lines
    parse_lines!(ex, p, x-> !( x in [_PARAMETERS[1], _CONSTPARAMETERS[1]] ))
    fill_unspecified_with_defaults(p)

    abstract_type = prepare_abstract_type(
        p.extras[_LINEAR[1]] ? :LinearDiffusion : :DiffusionProcess,
        p.dims,
        p.extras[_ELTYPE[1]],
        p.extras[_STATESPACE[1]],
    )
    #struct_def, struct_body, struct_const = createstruct(abstract_type, p)
    struct_def = createstruct(abstract_type, p)
    #struct_def = Expr(
    #    :macrocall,
    #    :(Core.@doc),
    #    join(repr(struct_body), "\n\n\t", repr(struct_const)),
    #    struct_def,
    #)

    add_diff_function!(p.fns, p)
    add_parameter_names_function!(p.fns, p)
    add_parameter_names_function!(p.fns, p, :const_parameter_names, p.constp, identity)
    #add_parameter_names_function!(p.fns, p, :var_parameter_names, p.nonconstp, identity)
    #add_end_point_info_function!(p.fns, p)
    add_end_point_info_names_function!(p.fns, p)

    eval(struct_def)
    for fn in p.fns
        eval(fn)
    end
    println("A new struct `$curly_name` has been defined.")
    println("To learn more about how to define an instance of this struct")
    println("please type in `?$name` and hit `ENTER`.")

    Meta.parse("import DiffusionDefinition.$name")
end

"""
    process_name

Process the name of a struct by returning the pure name (without template
arguments), the name with template arguments and a list of template arguments.
"""
function process_name end
process_name(name::Symbol) = name, name, Any[]
process_name(name::Expr) = name.args[1], name, name.args[2:end]

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
    _symbol_in(line, _CONSTPARAMETERS) && return _CONSTPARAMETERS[1]
    _symbol_in(line, _ENDPOINTINFO) && return _ENDPOINTINFO[1]
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
function parse_line!(
        ::K, line, p
    ) where K <: Union{Val{:parameters}, Val{:const_parameters}}
    param_name = line.args[1]

    if typeof(param_name) <: Symbol
        parse_param_single_name(line, p, K==Val{:const_parameters})
    else
        @assert param_name.head == :tuple
        parse_param_multi_names(line, p, K==Val{:const_parameters})
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
function parse_param_single_name(line, p, const_param=false)
    name_stem, disambig_idx, generic_name = get_name_stem(
        line.args[1], p.parameters
    )

    data_type = line.args[2]

    if _is_datatype(data_type, p)
        num_params = 1
    else
        @assert data_type.head == :tuple
        num_params, data_type = data_type.args[1], data_type.args[2]
        generic_name = true
    end

    for i in 1:num_params
        name = ( generic_name ? Symbol(name_stem, disambig_idx+i) : name_stem )
        append!(p.parameters, [(name, data_type, name_stem, disambig_idx+i)])
        push!( (const_param ? p.constp : p.nonconstp), name)
    end
end

"""
    _is_datatype(sym, p)

Utility function that checks whether `sym` is a datatype. It returns true if
`sym` is either an in-built datatype (say float, StaticArray etc.) or if it is a
template argument.
"""
function _is_datatype(sym, p)
    (sym in p.template_args) && return true
    typeof(sym) <: Expr && sym.head == :tuple && return false
    typeof(eval(sym)) <: DataType
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
function parse_param_multi_names(line, p, const_param=false)
    names = line.args[1].args
    data_types = line.args[2]

    if _is_datatype(data_types, p)
        data_types = repeat([data_types], length(names))
    else
        @assert data_types.head == :tuple
        if _is_datatype(data_types.args[1], p)
            @assert length(data_types.args) == length(names)
            @assert all( map(x->_is_datatype(x, p), data_types.args) )
            data_types = data_types.args
        else
            @assert typeof(eval(data_types.args[1])) <: Number
            @assert length(data_types.args) == 2
            @assert eval(data_types.args[1]) == length(names)
            @assert _is_datatype(data_types.args[2], p)
            data_types = repeat([data_types.args[2]], length(names))
        end
    end

    current_param_names = [param[1] for param in p.parameters]
    for (n,d) in zip(names, data_types)
        @assert !(n in current_param_names)
        append!(p.parameters, [(n,d,n,1)])
        push!( (const_param ? p.constp : p.nonconstp), n)
    end
end

function add_parameter_names_function!(
        fns, p, f_name=:parameter_names, coll=p.parameters, λ=(x->x[1])
    )

    fn_def = Expr(
        :call,
        f_name,
        Expr(
            :(::),
            Expr(
                :curly,
                :Type,
                Expr(
                    :(<:),
                    p.name,
                )
            )
        )
    )
    _, used_names = organize_end_point_info(p)

    fn_body = Expr(
        :tuple,
        [QuoteNode(λ(x)) for x in coll]...,
        [QuoteNode(name) for name in used_names]...,
    )
    push!(fns, Expr(:(=), fn_def, fn_body))
end


#------------------------------------------------------------------------------#
#
#                      For parsing: END POINT INFORMATION
#
#------------------------------------------------------------------------------#
"""
    parse_line!(::Val{:aux_info}, line, p)

Parse a line that defines parameters of the diffusion. The line must be in a
format:
    name --> parameter-description
"""
function parse_line!(::Val{:aux_info}, line, p)
    name = line.args[1]
    for i in 1:6
        _symbol_in(name, _ENDPOINTEXTRA[i]) && (name = _ENDPOINTEXTRA[i][1])
    end

    p.extras[name] = line.args[2]
end


function add_end_point_info_names_function!(fns, p)
    fn_def = Expr(
        :call,
        :end_point_info_names,
        Expr(
            :(::),
            Expr(
                :curly,
                :Type,
                Expr(
                    :(<:),
                    p.name,
                )
            )
        )
    )

    _, used_names = organize_end_point_info(p)

    fn_body = Expr(
        :tuple,
        [QuoteNode(name) for name in used_names]...,
    )

    push!(fns, Expr(:(=), fn_def, fn_body))
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

cleanup_param_names!(num::Number, t::Symbol, x::Symbol, params) = num

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
    _symbol_in(name, _SPARSEDIFF) && (name = _SPARSEDIFF[1])
    _symbol_in(name, _DIAGONALDIFF) && (name = _DIAGONALDIFF[1])
    _symbol_in(name, _SPARSEBMAT) && (name = _SPARSEBMAT[1])
    _symbol_in(name, _DIAGONALBMAT) && (name = _DIAGONALBMAT[1])
    _symbol_in(name, _LINEAR) && (name = _LINEAR[1])
    _symbol_in(name, _ELTYPE) && (name = _ELTYPE[1])

    p.extras[name] = line.args[2]
end

"""
    add_diff_function!(fns, p)

Add a definition of a function `consdiff` that indicates if the diffusion
coefficient is constant
"""
function add_diff_function!(fns, p)
    to_add = [
        (:constdiff, _CONSTDIFF[1]),
        (:sparsediff, _SPARSEDIFF[1]),
        (:diagonaldiff, _DIAGONALDIFF[1]),
        (:sparseBmat, _SPARSEBMAT[1]),
        (:diagonalBmat, _DIAGONALBMAT[1])
    ]
    for (name, _NAME) in to_add
        fn_def = Expr(:call,
            name,
            Expr(:(::),
                :P,
                p.name
            )
        )
        push!(fns, Expr(:(=), fn_def, p.extras[_NAME]))
    end
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
    _symbol_in(dim_of_what, _WIENER) && (p.dims[_WIENER[1]] = eval(dim))
    _symbol_in(dim_of_what, _PROCESS) && (p.dims[_PROCESS[1]] = eval(dim))
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
        (_SPARSEDIFF[1], false),
        (_DIAGONALDIFF[1], false),
        (_SPARSEBMAT[1], false),
        (_DIAGONALBMAT[1], false),
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
    Expr(:curly,
        stem,
        data_type,
        dims[_PROCESS[1]],
        dims[_WIENER[1]],
        state_restr,
    )
end

"""
    createstruct(abstract_type, name, params)

Create code that defines a struct defining a diffusion process.
"""
function createstruct(abstract_type, p)
    param_vec = map(p.parameters) do (name, data_type, _, _)
        Expr(:(::), name, data_type)
    end
    end_point_info_vec, end_point_info_used_names = organize_end_point_info(p)

    _new = (
        p.template_args == Any[] ?
        :new :
        Expr(:curly, :new, p.template_args...,)
    )

    xT_default = (
        length(end_point_info_vec)>0 ?
        [
            Expr(:kw,
                end_point_info_vec[end],
                Expr(
                    :call,
                    :zero,
                    end_point_info_vec[end].args[2],
                    p.dims[_PROCESS[1]],
                    Expr(
                        :call,
                        Expr(
                            :.,
                            :DiffusionDefinition,
                            :(:ismutable),
                        ),
                        end_point_info_vec[end].args[2],
                    )
                ),
            )
        ] :
        []
    )


    constructor_def = add_where_decorator(
        Expr(
            :call,
            p.name,
            param_vec...,
            end_point_info_vec[1:end-1]...,
            xT_default...
        ),
        p
    )

    named_constructor_def = add_where_decorator(
        Expr(
            :call,
            p.name,
            Expr(
                :parameters,
                param_vec...,
                end_point_info_vec...,
            )
        ),
        p
    )

    partial_clone_constructor_def = add_where_decorator(
        Expr(
            :call,
            p.name,
            Expr(
                :parameters,
                param_vec...,
            ),
            end_point_info_vec[1:end-1]...,
            xT_default...
        ),
        p
    )

    build_costructor(x)= Expr(
        :function, # constructor
        x,
        Expr(:block,
            Expr(:call,
                _new,
                map(x->x[1], p.parameters)...,
                end_point_info_used_names...,
            ),
        ), # body of the constructor
    ) # constructor

    struct_def = Expr(:struct,
        true, # immutable
        Expr(:<:,
            p.curly_name,
            abstract_type,
        ), # first line defining the struct
        Expr(:block,
            param_vec..., # parameters
            end_point_info_vec..., # extra info about end-points
            build_costructor(constructor_def),
            build_costructor(named_constructor_def),
            build_costructor(partial_clone_constructor_def),
        ), # body defining the struct
    )
    struct_def#, Expr(:block, param_vec..., end_point_info_vec...,), constructor_def
end

function add_where_decorator(expr, p)
    (
        p.template_args == Any[] ?
        expr :
        Expr(:where,
            expr,
            p.template_args...,
        )
    )
end


function organize_end_point_info(p)
    end_point_info_names = [:t0, :T, :v0, :x0, :vT, :xT]
    end_point_info_vec = Any[]
    end_point_info_used_names = Symbol[]
    if any([haskey(p.extras, eio) for eio in end_point_info_names])
        for name in [:t0, :T, :v0, :x0, :vT]
            haskey(p.extras, name) && (
                push!(end_point_info_vec, Expr(:(::), name, p.extras[name]));
                push!(end_point_info_used_names, name)
            )
        end
        # the space for the exact end-point must always be there due to blocking
        default_datatype = find_state_datatype(
            p.extras[_ELTYPE[1]], # eltype
            p.dims[_PROCESS[1]], # dim of the process
        )
        push!(
            end_point_info_vec,
            Expr(:(::), :xT, get(p.extras, :xT, default_datatype))
        )
        push!(end_point_info_used_names, :xT)
    end
    end_point_info_vec, end_point_info_used_names
end


function find_state_datatype(_eltype, dim_process)
    dim_process < 9 && return SVector{dim_process,eval(_eltype)}
    return Vector{eval(_eltype)}
end
