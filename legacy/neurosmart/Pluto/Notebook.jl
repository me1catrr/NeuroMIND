### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 60620496-6fc4-45df-a2fc-eca64e76a52f
begin
	using PyMNE
	using Plots
	using CSV
	using DataFrames
	using Serialization
	using Dates 
	using DSP
	using Serialization
	using Statistics
	using StatsBase
	using StatsPlots
	using FFTW
	using UnfoldMakie;
	using CairoMakie
	using DataFramesMeta
	using UnfoldSim
	using Unfold
	using MakieThemes
	using TopoPlots;
end

# ╔═╡ f36b8ac2-d883-4d32-a7dc-2072e987165d
#=
begin
# 1) Montaje base
montage = PyMNE.channels.make_standard_montage("standard_1020")

# 2) Labels y posiciones 3D
labels_all = pyconvert(Vector{String}, montage.ch_names)
ch_pos = pyconvert(Dict{String,Any}, montage.get_positions()["ch_pos"])

# 3) Elige tus 32 canales
wanted = [
    "Fp1","Fp2",
    "F7","F3","Fz","F4","F8",
    "FC5","FC1","FC2","FC6",
    "T7","C3","Cz","C4","T8",
    "CP5","CP1","CP2","CP6",
    "P7","P3","Pz","P4","P8",
    "PO9","O1","Oz","O2","PO10","TP9","TP10"
]

labels = [l for l in labels_all if l in wanted]

# 4) Posiciones 3D -> 2D
pos3d = hcat([ch_pos[l] for l in labels]...)
positions = to_positions(pos3d)

# 5) Tus datos: un valor por canal
# reemplaza esto con tus betas/ERP/voltajes/etc.
values = randn(length(labels))

# 6) Topoplot estilo UnfoldMakie
plot_topoplot(
    values;
    labels = labels,
    positions = Point2f.(positions),
    visual = (; label_text = true, label_scatter = false),
    axis = (; xlabel = "")
)
end
=#

# ╔═╡ a73289ad-5d2c-4452-b19e-8a0ef3e42836
#=
begin
	x = 0:0.01:10
	y = sin.(x)

	fig = Figure()
	ax = Axis(fig[1,1], xlabel="x", ylabel="sin(x)")

	lines!(ax, x, y)

	fig
end
=#