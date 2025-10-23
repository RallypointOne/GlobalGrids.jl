### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# ╔═╡ beddf43c-afef-11f0-82e4-cff7e7077134
begin 
	using Pkg
	Pkg.activate(@__DIR__)
	using DGGS, BenchmarkTools, DimensionalData, YAXArrays, PlutoUI, GLMakie

	PlutoUI.TableOfContents()
end

# ╔═╡ 08bea61d-f6ab-40c4-887b-122d54fcf6dd
md"# Benchmark to_cell and to_geo"

# ╔═╡ f8216103-5901-408f-8cb5-303521994ca8
@benchmark to_cell(120.8, -80, 8)

# ╔═╡ 3b455d44-39ef-491d-8c8d-3b35694e838f
@benchmark to_geo($(to_cell(120.8, -80, 8)))

# ╔═╡ 90675a51-ddbd-4829-b3ad-f55273a59a90
md"# Plot"

# ╔═╡ c0a26896-02a2-4b5b-8d23-fc5b796a5cb7
begin 
	lon_range = X(180:-1:-180)
	lat_range = Y(90:-1:-90)
	geo_data = [exp(cosd(lon)) + 3(lat / 90) for lon in lon_range, lat in lat_range]
	geo_array = YAXArray((lon_range, lat_range), geo_data)
	plot(geo_array)
end

# ╔═╡ 57035de2-2770-4fb9-b9c5-6647befe0bf0
begin
	resolution = 4
	dggs_pyramid = to_dggs_pyramid(geo_array, resolution, "EPSG:4326")
end

# ╔═╡ a51f47be-5f1b-4ee1-9ff2-9bdc2892121c
dggs_array = dggs_pyramid[2].layer1

# ╔═╡ 4add4abd-0d6e-4ce9-9c15-d6389d35e47b
plot(dggs_array)

# ╔═╡ Cell order:
# ╠═beddf43c-afef-11f0-82e4-cff7e7077134
# ╠═08bea61d-f6ab-40c4-887b-122d54fcf6dd
# ╠═f8216103-5901-408f-8cb5-303521994ca8
# ╠═3b455d44-39ef-491d-8c8d-3b35694e838f
# ╟─90675a51-ddbd-4829-b3ad-f55273a59a90
# ╠═c0a26896-02a2-4b5b-8d23-fc5b796a5cb7
# ╠═57035de2-2770-4fb9-b9c5-6647befe0bf0
# ╠═a51f47be-5f1b-4ee1-9ff2-9bdc2892121c
# ╠═4add4abd-0d6e-4ce9-9c15-d6389d35e47b
