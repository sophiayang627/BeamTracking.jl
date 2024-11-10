module Misc
using ..GTPSA: @FastGTPSA!, GTPSA, evaluate!, TPS
import ..BeamTracking: track!
using ..BeamTracking
export track!

Base.@kwdef struct Taylor{T<:TPS,U<:Union{Nothing,Quaternion{<:T}}}
  v::Coord{T}
  q::U
end

function track!(beam::Beam, ele::Misc.Taylor)

end


end