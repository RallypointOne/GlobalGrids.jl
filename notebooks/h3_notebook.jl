### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ b0305888-b8c5-11f0-b4a3-439b3025b969
begin 
	using Pkg 
	Pkg.activate(joinpath(@__DIR__, ".."))
	push!(LOAD_PATH, joinpath(@__DIR__, ".."))
	unique!(LOAD_PATH)

	using GlobalGrids, GLMakie, GLMakie 
	import GeoInterface as GI
	import GlobalGrids as GG
end

# ╔═╡ f7e66f01-a090-4ce8-9902-4017dd5777da


# ╔═╡ Cell order:
# ╠═b0305888-b8c5-11f0-b4a3-439b3025b969
# ╠═f7e66f01-a090-4ce8-9902-4017dd5777da
