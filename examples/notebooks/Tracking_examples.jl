### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# ╔═╡ 2a01b3c0-1bcb-11f0-29c1-35727198cb30
begin
	# Activate beam tracking package
	using Pkg;
	cd(@__DIR__)
	cd("../..")
	Pkg.activate(".")
	Pkg.add("LaTeXStrings")

end

# ╔═╡ e4c10317-3276-4a7a-a270-06c03ab6edd5
using BeamTracking, Plots, LaTeXStrings# , Unitful, PhysicalConstants

# ╔═╡ 260c8d25-df28-4f79-8af5-49cdaba00b18
md"# Tracking examples"

# ╔═╡ ca4b6a03-5862-4add-926e-72dc6c9715be
## Press an "eye" symbol on the left to show/hide cell code

# ╔═╡ 1db55000-c591-4ffb-a039-490da5ca1a43
md"This is also a way to input code blocks into the **Markdown** cells"

# ╔═╡ bf3c07f5-bba3-4355-b51e-62b695a840d1
md"In future we want to just have 
```Julia
using BeamTracking 
```
and Pluto downloads BeamTracking and prerequisites 
"

# ╔═╡ 44db8dc5-4bd8-4e09-80ff-cf0843beda99
#import  PhysicalConstants.CODATA2022: c_0 as c, m_e as m

# ╔═╡ ade360c6-af69-4df1-aa29-3d83a64a7d18
md"Define some bunch parameters"

# ╔═╡ d0934715-3371-47ae-bf14-59bebc79ba15
begin
	e_minus = Species("electron")
	mec2 = massof(e_minus)
	#Kinetic energy
	ek = 5e3
	# Beta - gamma
	bg = sqrt(ek / mec2 * (ek / mec2 + 2))
end

# ╔═╡ 08ec929a-ae56-44c1-b538-a4502389535e
md"Create bunch"

# ╔═╡ 8b51f3f5-0a53-4dd5-8cb9-5a827731cdcf
begin
	NParticles = 256
	xi  = -2e-3 .+ 4e-3 	* randn(NParticles)
	pxi = -7.5e-4 .+ 1.5e-3 * randn(NParticles)
	yi  = -1e-3 .+ 2e-3 		* randn(NParticles)
	pyi = -3e-4 .+ 6e-4 		* randn(NParticles)
	zi  = zeros(NParticles)
	pzi = -1e-3 .+ 2e-3*randn(NParticles)

	bunch = Bunch(
		species = e_minus, 
		beta_gamma_ref = bg,       
		x = copy(xi), px = copy(pxi), 
		y = copy(yi), py = copy(pyi), 
		z = copy(zi), pz = copy(pzi)
	)
	plot(
		scatter(bunch.v.x, bunch.v.px, xlabel = L"x", ylabel = L"p_x"),
		scatter(bunch.v.y, bunch.v.py, xlabel = L"x", ylabel = L"p_x"),
		scatter(bunch.v.z, bunch.v.pz, xlabel = L"z", ylabel = L"p_z"),
		layout = (3,1), size=(500,900)
	)

end

# ╔═╡ ea1763aa-b397-4cd9-a751-d3fd5fe08bb7
md"Create some drifts"

# ╔═╡ b047b953-d21f-43a8-a619-68c9256d517d
begin
	ld1 = 0.15;  # m
	ld2 = 0.75;  # m
	ld3 = 2.00;  # m
	dr1 = MatrixKick.Drift(L = ld1)
	dr2 = MatrixKick.Drift(L = ld2)
	dr3 = MatrixKick.Drift(L = ld3)
end

# ╔═╡ 116876c3-8daa-43d3-92f2-069eea14c186
md"Track and plot"

# ╔═╡ 96d8a130-d3b3-421b-99c0-736302f85a22
begin
	track!(bunch, dr1);
	plot(
		scatter(bunch.v.x, bunch.v.px, xlabel = L"x", ylabel = L"p_x"),
		scatter(bunch.v.y, bunch.v.py, xlabel = L"x", ylabel = L"p_x"),
		scatter(bunch.v.z, bunch.v.pz, xlabel = L"z", ylabel = L"p_z"),
		layout = (3,1), size=(500,900)
	)
end

# ╔═╡ Cell order:
# ╟─260c8d25-df28-4f79-8af5-49cdaba00b18
# ╠═ca4b6a03-5862-4add-926e-72dc6c9715be
# ╟─1db55000-c591-4ffb-a039-490da5ca1a43
# ╠═2a01b3c0-1bcb-11f0-29c1-35727198cb30
# ╟─bf3c07f5-bba3-4355-b51e-62b695a840d1
# ╠═e4c10317-3276-4a7a-a270-06c03ab6edd5
# ╠═44db8dc5-4bd8-4e09-80ff-cf0843beda99
# ╟─ade360c6-af69-4df1-aa29-3d83a64a7d18
# ╠═d0934715-3371-47ae-bf14-59bebc79ba15
# ╟─08ec929a-ae56-44c1-b538-a4502389535e
# ╠═8b51f3f5-0a53-4dd5-8cb9-5a827731cdcf
# ╟─ea1763aa-b397-4cd9-a751-d3fd5fe08bb7
# ╠═b047b953-d21f-43a8-a619-68c9256d517d
# ╠═116876c3-8daa-43d3-92f2-069eea14c186
# ╠═96d8a130-d3b3-421b-99c0-736302f85a22
