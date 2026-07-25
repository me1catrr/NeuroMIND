#!/usr/bin/env julia
# Regenera figuras de publicación EC desde CSV ya existentes (sin re-análisis).
# Uso: julia --project=. scripts/regenerate_group_figures.jl [EC|EO|both]

include(joinpath(@__DIR__, "..", "src", "viz", "GroupVizCommon.jl"))
using .GroupVizCommon
using CSV, DataFrames, Statistics

const PROJ = dirname(@__DIR__)
const ROOT = resolve_results_root(PROJ)

function regen_long(cond::String)
    dir = joinpath(ROOT, "longitudinal", cond)
    isdir(dir) || return @warn "No existe $dir"
    figs = joinpath(dir, "figures"); mkpath(figs)
    println("── Longitudinal $cond → $figs")
    sm = safe_csv(joinpath(dir, "subject_band_means.csv"))
    nrow(sm) > 0 && save_paired_means(joinpath(figs, "paired_mean_wpli_by_band.png"), sm)
    for b in BAND_ORDER
        r1 = read_mat_csv(joinpath(dir, "longitudinal_connectivity_t1_$(b).csv"))
        r2 = read_mat_csv(joinpath(dir, "longitudinal_connectivity_t2_$(b).csv"))
        rd = read_mat_csv(joinpath(dir, "longitudinal_difference_$(b).csv"))
        (r1 === nothing || r2 === nothing) && continue
        ch, Wt1, Wt2 = r1[1], r1[2], r2[2]
        Wd = rd === nothing ? (Wt2 .- Wt1) : rd[2]
        lim_ab = max(maximum(Wt1), maximum(Wt2), 1e-12)
        save_heatmap(joinpath(figs, "heatmap_t1_$(b).png"), Wt1, ch, "wPLI T1 — $b ($cond)";
                     colorrange=(0.0, lim_ab))
        save_heatmap(joinpath(figs, "heatmap_t2_$(b).png"), Wt2, ch, "wPLI T2 — $b ($cond)";
                     colorrange=(0.0, lim_ab))
        save_heatmap(joinpath(figs, "heatmap_delta_$(b).png"), Wd, ch, "Δ wPLI (T2−T1) — $b ($cond)";
                     diverging=true, colorbar_label="Δ wPLI")
        save_heatmap_triplet(joinpath(figs, "heatmap_triplet_$(b).png"), Wt1, Wt2, Wd, ch,
            ("T1 — $b ($cond)", "T2 — $b ($cond)", "Δ (T2−T1) — $b ($cond)"))
        sdf = safe_csv(joinpath(dir, "longitudinal_statistics_$(b).csv"))
        if nrow(sdf) > 0
            for c in (:ch_a, :ch_b); hasproperty(sdf, c) && (sdf[!, c] = String.(sdf[!, c])); end
            sig = filter(r -> Float64(r.q_value) < 0.05, sdf)
            if nrow(sig) > 0
                save_topo_network(joinpath(figs, "sig_network_$(b).png"), ch,
                    [NamedTuple(r) for r in eachrow(sig)], "Edges FDR — $b ($cond)")
            else
                top = first(sort(sdf, :effect_d; by=abs, rev=true), min(20, nrow(sdf)))
                save_topo_network(joinpath(figs, "explore_network_topN_$(b).png"), ch,
                    [NamedTuple(r) for r in eachrow(top)],
                    "Top-20 |dz| (exploratorio) — $b ($cond)")
            end
        end
        println("  $b OK ($(length(ch)) ch)")
    end
    bp = safe_csv(joinpath(dir, "tables", "spectral", "band_power_delta_statistics.csv"))
    if nrow(bp) > 0
        for b in BAND_ORDER
            sub = filter(r -> string(r.band) == b, bp)
            nrow(sub) == 0 && continue
            chs = String.(sub.channel)
            dvs = Float64.(sub.diff)
            save_topo_delta(joinpath(figs, "topo_delta_bandpower_$(b).png"), chs, dvs,
                            Dict{String,Tuple{Float64,Float64}}(),
                            "Δ band power — $b ($cond)"; colorbar_label="Δ power (T2−T1)")
        end
    end
end

function regen_trans(cond::String)
    dir = joinpath(ROOT, "transversal", cond)
    isdir(dir) || return @warn "No existe $dir"
    figs = joinpath(dir, "figures"); mkpath(figs)
    println("── Transversal $cond → $figs")
    sm = safe_csv(joinpath(dir, "subject_band_means.csv"))
    nrow(sm) > 0 && save_group_means(joinpath(figs, "group_mean_wpli_by_band.png"), sm)
    for b in BAND_ORDER
        rc = read_mat_csv(joinpath(dir, "group_connectivity_control_$(b).csv"))
        rm = read_mat_csv(joinpath(dir, "group_connectivity_ms_$(b).csv"))
        rd = read_mat_csv(joinpath(dir, "group_difference_$(b).csv"))
        (rc === nothing || rm === nothing) && continue
        ch, Wctrl, Wms = rc[1], rc[2], rm[2]
        Wd = rd === nothing ? (Wms .- Wctrl) : rd[2]
        lim_ab = max(maximum(Wctrl), maximum(Wms), 1e-12)
        save_heatmap(joinpath(figs, "heatmap_control_$(b).png"), Wctrl, ch, "wPLI Control — $b ($cond)";
                     colorrange=(0.0, lim_ab))
        save_heatmap(joinpath(figs, "heatmap_ms_$(b).png"), Wms, ch, "wPLI MS — $b ($cond)";
                     colorrange=(0.0, lim_ab))
        save_heatmap(joinpath(figs, "heatmap_diff_$(b).png"), Wd, ch, "Δ wPLI (MS−Control) — $b ($cond)";
                     diverging=true, colorbar_label="Δ wPLI")
        save_heatmap_triplet(joinpath(figs, "heatmap_triplet_$(b).png"), Wctrl, Wms, Wd, ch,
            ("Control — $b ($cond)", "MS — $b ($cond)", "Δ (MS−Ctrl) — $b ($cond)"))
        sdf = safe_csv(joinpath(dir, "group_statistics_$(b).csv"))
        if nrow(sdf) > 0
            for c in (:ch_a, :ch_b); hasproperty(sdf, c) && (sdf[!, c] = String.(sdf[!, c])); end
            sig = filter(r -> Float64(r.q_value) < 0.05, sdf)
            if nrow(sig) > 0
                save_topo_network(joinpath(figs, "sig_network_$(b).png"), ch,
                    [NamedTuple(r) for r in eachrow(sig)], "Edges FDR — $b ($cond)")
            else
                top = first(sort(sdf, :effect_d; by=abs, rev=true), min(20, nrow(sdf)))
                save_topo_network(joinpath(figs, "explore_network_topN_$(b).png"), ch,
                    [NamedTuple(r) for r in eachrow(top)],
                    "Top-20 |d| (exploratorio) — $b ($cond)")
            end
        end
        println("  $b OK ($(length(ch)) ch)")
    end
    bp = safe_csv(joinpath(dir, "tables", "spectral", "band_power_group_statistics.csv"))
    if nrow(bp) > 0
        for b in BAND_ORDER
            sub = filter(r -> string(r.band) == b, bp)
            nrow(sub) == 0 && continue
            chs = String.(sub.channel)
            dvs = Float64.(sub.diff)
            save_topo_delta(joinpath(figs, "topo_diff_bandpower_$(b).png"), chs, dvs,
                            Dict{String,Tuple{Float64,Float64}}(),
                            "Δ band power (MS−Control) — $b ($cond)";
                            colorbar_label="Δ power (MS−Control)")
        end
    end
end

conds = ["EC"]
if length(ARGS) > 0
    a = uppercase(ARGS[1])
    conds = a == "BOTH" ? ["EC", "EO"] : [a]
end

for c in conds
    regen_long(c)
    regen_trans(c)
end
println("✓ Figuras regeneradas.")
