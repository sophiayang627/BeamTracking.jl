#=

Temp code for AAPC

=#

Base.@kwdef struct Species
  name::String
end
massof(::Species) = 510998.95069*u"eV/c^2"
chargeof(::Species) = -1*u"q"

c_light = 2.99792458e8*u"m/s"