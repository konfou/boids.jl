## HIDDEN SECTION
#

using LinearAlgebra
using Random
using Makie, GLMakie, AbstractPlotting
#import Plots
AbstractPlotting.inline!(false)

## VISIBLE SECTION
#

mutable struct Boid{T <: Float64}
  id::Int; pos::Array{T,1}; vel::Array{T,1}
  dv::T; av::T; ds::T
  w0::T; ws::T; wa::T; wc::T
  vmax::T; vmin::T
  bmin::Array{T,1}; bmax::Array{T,1}
  ww::T; wb::T; wp::T; we::T
  nc::Int; dp::T; de::T;
end

dis(bi, bj) = norm(bi.pos - bj.pos)
ang(bi, bj) = acos(dot(bi.vel, bj.pos)/(norm(bi.vel)*norm(bj.pos)))
dim(b) = length(b.pos)

## HIDDEN SECTION
#

macro reload()
    :( include("boids.jl"); println("Reloaded.") )
end

mutable struct World
  xy; sc; ms
  World() = new()
end

include("models/basic.jl")
include("models/modded.jl")
include("models/pureva.jl")
