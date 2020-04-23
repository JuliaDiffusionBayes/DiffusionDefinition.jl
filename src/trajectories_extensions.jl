

Trajectories.trajectory(tt, v::Type, D, ::Val{true}) = trajectory(tt, v, D)
Trajectories.trajectory(tt, v::Type, D, ::Val{false}) = trajectory(tt, v)

"""
    Trajectories.trajectory(tt, v::Type, D::Number)

Create a Trajectory with mutable states of dimension `D` along the time
collection `tt`.
"""
function Trajectories.trajectory(tt, v::Type, D::Number)
    trajectory(collect(tt), [zeros(eltype(v),D) for _ in tt])
end

"""
    Trajectories.trajectory(tt, v::Type)

Create a Trajectory with immutable states along the time collection `tt`.
"""
function Trajectories.trajectory(tt, v::Type)
    ismutable(v)::Val{false}
    trajectory(collect(tt), zeros(v, length(tt)))
end


function Trajectories.trajectory(
        tt::AbstractArray{<:AbstractArray},
        P::DiffusionProcess,
        v::Type=default_type(P),
        w::Type=default_wiener_type(P),
    )
    (
        process = [_process_traj(t, P, v) for t in tt],
        wiener = [_wiener_traj(t, P, w) for t in tt],
    )
end

function Trajectories.trajectory(
        tt,
        P::DiffusionProcess,
        v::Type=default_type(P),
        w::Type=default_wiener_type(P)
    )
    (
        process = _process_traj(tt, P, v),
        wiener = _wiener_traj(tt, P, w)
    )
end

function _process_traj(tt, ::DiffusionProcess{T,DP}, v::Type) where {T,DP}
    trajectory(tt, v, eval(DP))
end

function _wiener_traj(tt, ::DiffusionProcess{T,DP,DW}, v::Type) where {T,DP,DW}
    trajectory(tt, v, eval(DW))
end


@recipe function f(tr::Trajectory, ::Val{:vs_time}; coords=1:length(tr.x[1]))
    for i in coords
        @series begin
            seriestype := :path
            tr.t, map(x->x[i], tr.x)
        end
    end
end

@recipe function f(tr::Trajectory, ::Val{:x_vs_y}; coords=1:length(tr.x[1]))
    @assert length(coords) == 2
    @series begin
        seriestype := :path
        map(x->x[coords[1]], tr.x), map(x->x[coords[2]], tr.x)
    end
end


@recipe function f(tr::Trajectory; vs_time=true, coords=1:length(tr.x[1]))
    plot(tr, vs_time ? Val(:vs_time) : ValVal(:x_vs_y); coords=coords)
end
