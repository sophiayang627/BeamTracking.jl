#=

Temp code for AAPC

=#

Base.@kwdef struct Species
  name::String
end
massof(::Species) = 510998.95069
chargeof(::Species) = -1

const c_light = 2.99792458e8