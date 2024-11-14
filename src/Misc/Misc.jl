module Misc
using ..GTPSA: @FastGTPSA!, GTPSA, evaluate!, TPS
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

Base.@kwdef struct Taylor{T<:TPS,U<:Union{Nothing,Quaternion{<:T}}}
  v::Coord{T}
  q::U
end

function track!(bunch::Bunch, ele::Misc.Taylor)

end


end