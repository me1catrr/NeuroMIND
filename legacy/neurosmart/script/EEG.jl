#!/usr/bin/env julia
# -*- coding: utf-8 -*-

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

# Configuración (mantengo el módulo actual)
include(joinpath(@__DIR__, "..", "config", "default_config.jl"))
using .DefaultConfig
cfg = DEFAULT_CONFIG

# Cargar módulo del proyecto
include(joinpath(@__DIR__, "..", "src", "modules", "EEG_Julia.jl"))
using .EEG_Julia

# Pipeline secuencial (flujo visible)
# NOTA: cada `run_*` debe poder ejecutarse también por separado (REPL/Pluto).
run_io(cfg)
run_filtering(cfg)
run_ica(cfg)
run_ica_cleaning(cfg)
run_segmentation(cfg)
run_baseline(cfg)
run_artifact_rejection(cfg)
run_baseline_2st(cfg)
run_fft(cfg)
run_csd(cfg)
run_wpli(cfg)