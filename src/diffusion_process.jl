using MacroTools, StaticArrays
import Base: lowercase
import Base: eltype

abstract type DiffusionProcess{T,DP,DW} end
dimension(d::DiffusionProcess{T,DP,DW}) where {T,DP,DW} = (
    process = DP,
    wiener = DW
)
eltype(d::DiffusionProcess{T}) where T = T

abstract type LinearDiffusion{T,DP,DW} <: DiffusionProcess{T,DP,DW} end

# Temporary definition
struct UnboundedDomain end

macro diffusion_process(name, ex::Expr, p...)
    parse_process(name, MacroTools.striplines(ex), p)
end

_DIMENSIONS = [:dimensions, :dim, :dims, :dimension]
_PARAMETERS = [:parameters, :param, :params]
_CONJUGATE = [:conjugate]
_EXTRA = [:additional, :extra]
_WIENER = [:wiener, :noise, :gaus, :gaussian]
_PROCESS = [:proc, :process, :state, :statespace]
_DOMAIN = [:domain, :domains, :statespace, :state]
_CONSTDIFF = [:constdiff, :constvola, :constdiffusivity, :constvolatility, :constantdiff, :constantvola]
_LINEAR = [:linear, :lineardiffusion]
_ELTYPE = [:eltype]
_NUMNONHYPO = [:num_non_hypo, :numnonhypo]
_PHI = [:phi, :ϕ, :φ]

lowercase(s::Symbol) = Symbol(lowercase(string(s)))
_symbol_in(s::Symbol, symbols) = lowercase(s) in symbols
_symbol_in(s::QuoteNode, symbols) = _symbol_in(eval(s), symbols)
_symbol_in(s::Any, ::Any) = false

_symbol_in(:dimensions, _DIMENSIONS)

function _update_label(line, current_label)
    _symbol_in(line, _DIMENSIONS) && return _DIMENSIONS[1]
    _symbol_in(line, _PARAMETERS) && return _PARAMETERS[1]
    _symbol_in(line, _CONJUGATE) && return _CONJUGATE[1]
    _symbol_in(line, _EXTRA) && return _EXTRA[1]
    current_label
end

function parse_line!(::Val{:parameters}, line, p)
    data_type = line.args[1]

    if typeof(data_type) <: Symbol
        parse_param_single_name(line, p)
    else
        @assert data_type.head == :tuple
        parse_param_multi_names(line, p)
    end
end

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
    end
end

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

function add_phi_function!(expr::Expr, fns::Vector{Expr}, t::Symbol, x::Symbol, p)
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

function parse_line!(::Val{:additional}, line, p)
    name = line.args[1]
    _symbol_in(name, _DOMAIN) && (name = _DOMAIN[1])
    _symbol_in(name, _CONSTDIFF) && (name = _CONSTDIFF[1])
    _symbol_in(name, _LINEAR) && (name = _LINEAR[1])
    _symbol_in(name, _ELTYPE) && (name = _ELTYPE[1])

    p.extras[name] = line.args[2]
end

function parse_line!(::Val{:dimensions}, line, p)
    dim_of_what, dim = line.args[1], line.args[2]
    _symbol_in(dim_of_what, _WIENER) && (p.dims[_WIENER[1]] = dim)
    _symbol_in(dim_of_what, _PROCESS) && (p.dims[_PROCESS[1]] = dim)
end


function parse_process(name , ex::Expr, ::Any)
    p = (
        name = name,
        parameters = Vector{Tuple{Symbol,Union{Symbol,Expr},Symbol,Int64}}(undef, 0),
        dims = Dict{Symbol,Int64}(),
        extras = Dict{Symbol, Any}(),
        fns = Vector{Expr}(undef, 0),
    )
    parse_lines!(ex, p, x->(x == _PARAMETERS[1]))
    parse_lines!(ex, p, x->(x != _PARAMETERS[1]))
    fill_unspecified_with_defaults(p)

    abstract_type = prepare_abstract_type(
        p.extras[_LINEAR[1]] ? :LinearDiffusion : :DiffusionProcess,
        p.dims,
        p.extras[_ELTYPE[1]]
    )
    struct_def = createstruct(abstract_type, name, p.parameters)
    eval(struct_def)
    for fn in p.fns
        eval(fn)
    end
    Meta.parse("import DiffusionDefinition.$name")
end

function parse_lines!(ex::Expr, p, condition)
    current_label = nothing
    for line in ex.args
        current_label = _update_label(line, current_label)
        typeof(line) <: Union{Symbol, QuoteNode} && (continue)
        @assert line.head == :-->
        condition(current_label) && parse_line!(Val{current_label}(), line, p)
    end
end

function fill_unspecified_with_defaults(p)
    fill_unspecified_with_defaults(Val{:dimensions}(), p)
    fill_unspecified_with_defaults(Val{:additional}(), p)
end

function fill_unspecified_with_defaults(::Val{:dimensions}, p)
    for key in [_WIENER[1], _PROCESS[1]]
        !haskey(p.dims, key) && (p.dims[key] = 1)
    end
end

function fill_unspecified_with_defaults(::Val{:additional}, p)
    defaults = [
        (_DOMAIN[1], UnboundedDomain()),
        (_CONSTDIFF[1], false),
        (_LINEAR[1], false),
        (_ELTYPE[1], Float64)
    ]
    for (key, default_val) in defaults
        !haskey(p.extras, key) && (p.extras[key] = default_val)
    end
end

function get_name_stem(name_stem::Symbol, parameters)
    if name_stem == :_
        name_stem = :p
    end
    disambig_idx = highest_idx_used(name_stem, parameters)
    name_stem, disambig_idx, (name_stem == :p)
end

function highest_idx_used(name_stem, params)
    idxs = [p[4] for p in params if name_stem == p[3]]
    length(idxs) == 0 && return 0
    maximum(idxs)
end

function prepare_abstract_type(stem, dims, data_type)
    wiener_dim, proc_dim = dims[_WIENER[1]], dims[_PROCESS[1]]
    "$stem{$data_type,$proc_dim,$wiener_dim}"
end

function createstruct(abstract_type, name, params)
    header = "struct $name <: $abstract_type\n"
    params_vec = map(params) do (name, data_type, _, _)
        "\t$name::$data_type\n"
    end
    Meta.parse(string(header, join(params_vec), "end\n"))
end
