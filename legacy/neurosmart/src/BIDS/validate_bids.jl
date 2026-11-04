#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/BIDS/validate_bids.jl
#
# VALIDACIÓN DE ESTRUCTURA BIDS
# ==============================
# Esta rutina valida de forma básica la estructura de un dataset EEG en formato
# BIDS, verificando tanto metadatos globales como la organización por sujeto/sesión.
#
# OBJETIVO:
#   - Detectar errores críticos que impiden el uso del dataset.
#   - Reportar advertencias que no bloquean, pero indican inconsistencias.
#
# ALCANCE:
#   - No reemplaza el validador oficial de BIDS.
#   - Cubre verificaciones frecuentes para pipelines EEG locales.
#
# FLUJO GENERAL:
#   1. Verifica la presencia de archivos obligatorios en la raíz
#   2. Valida dataset_description.json
#   3. Valida participants.tsv
#   4. Valida estructura de participantes
#   5. Valida archivos de tareas
#
# NOTA: Esta rutina es una herramienta de validación básica y no cubre todas
# las validaciones posibles de un dataset BIDS completo.

# Listado de librerías
using JSON              # Lectura y escritura de archivos JSON
using DataFrames        # Manipulación de datos tabulares
using CSV               # Lectura y escritura de archivos CSV
using Glob              # Buscar archivos en directorios

# ------------------------------------------------------------------------------------
# 1. VALIDACIÓN DE ESTRUCTURA BIDS
# ------------------------------------------------------------------------------------  
# Función principal para validar la estructura BIDS
# Parámetros:
#   - bids_dir: directorio raíz del dataset BIDS
# Retorna:
#   - errors: lista de errores críticos encontrados
#   - warnings: lista de advertencias encontradas
function validate_bids_structure(bids_dir)
    println("🔍 Validando estructura BIDS...")

    # Archivos mínimos esperados en la raíz de un dataset BIDS.
    required_files = ["dataset_description.json", "participants.tsv", "participants.json", "README"]
    for file in required_files
        if !isfile(joinpath(bids_dir, file))
            push!(errors, "❌ Archivo obligatorio faltante: $file")
        else
            println("✅ $file encontrado")
        end
    end
    
    # dataset_description.json define identidad del dataset y versión de BIDS.
    if isfile(joinpath(bids_dir, "dataset_description.json"))
        try
            dataset_desc = JSON.parsefile(joinpath(bids_dir, "dataset_description.json"))
            required_fields = ["Name", "BIDSVersion", "DatasetType"]
            for field in required_fields
                if !haskey(dataset_desc, field)
                    push!(errors, "❌ Campo obligatorio faltante en dataset_description.json: $field")
                end
            end
            println("✅ dataset_description.json válido")
        catch e
            push!(errors, "❌ Error parseando dataset_description.json: $e")
        end
    end
    
    # participants.tsv debe contener, al menos, la columna participant_id.
    if isfile(joinpath(bids_dir, "participants.tsv"))
        try
            participants = CSV.read(joinpath(bids_dir, "participants.tsv"), DataFrame)
            if !("participant_id" in names(participants))
                push!(errors, "❌ Columna 'participant_id' faltante en participants.tsv")
            else
                println("✅ participants.tsv válido ($(nrow(participants)) participantes)")
            end
        catch e
            push!(errors, "❌ Error leyendo participants.tsv: $e")
        end
    end
    
    # Se consideran participantes las carpetas con prefijo sub-.
    participant_dirs = filter(x -> startswith(x, "sub-"), readdir(bids_dir))
    println("📁 Encontrados $(length(participant_dirs)) directorios de participantes")
    
    for sub_dir in participant_dirs
        sub_path = joinpath(bids_dir, sub_dir)
        
        # Verificar que sea un directorio
        if !isdir(sub_path)
            push!(warnings, "⚠️ $sub_dir no es un directorio")
            continue
        end
        
        # Dentro de cada participante se esperan sesiones con prefijo ses-.
        session_dirs = filter(x -> startswith(x, "ses-"), readdir(sub_path))
        if isempty(session_dirs)
            push!(warnings, "⚠️ $sub_dir no tiene sesiones")
            continue
        end
        
        for ses_dir in session_dirs
            ses_path = joinpath(sub_path, ses_dir)
            eeg_path = joinpath(ses_path, "eeg")
            
            if !isdir(eeg_path)
                push!(warnings, "⚠️ $sub_dir/$ses_dir no tiene directorio eeg/")
                continue
            end
            
            # En cada sesión EEG, se valida la presencia de metadatos y señal.
            eeg_files = readdir(eeg_path)
            eeg_json_files = filter(x -> endswith(x, "_eeg.json"), eeg_files)
            eeg_data_files = filter(x -> endswith(x, "_eeg.eeg") || endswith(x, "_eeg.vhdr"), eeg_files)
            
            if isempty(eeg_json_files)
                push!(warnings, "⚠️ $sub_dir/$ses_dir/eeg/ no tiene archivos _eeg.json")
            end
            
            if isempty(eeg_data_files)
                push!(warnings, "⚠️ $sub_dir/$ses_dir/eeg/ no tiene archivos de datos EEG")
            end
            
            # Cada sidecar JSON se valida de forma independiente para aislar errores.
            for json_file in eeg_json_files
                json_path = joinpath(eeg_path, json_file)
                try
                    eeg_json = JSON.parsefile(json_path)
                    
                    # Campos mínimos útiles para procesamiento automático EEG.
                    required_eeg_fields = ["TaskName", "SamplingFrequency", "EEGChannelCount"]
                    for field in required_eeg_fields
                        if !haskey(eeg_json, field)
                            push!(warnings, "⚠️ $json_file falta campo: $field")
                        end
                    end
                    
                    # Consistencia semántica: TaskName del JSON vs nombre del archivo.
                    if haskey(eeg_json, "TaskName")
                        task_name = eeg_json["TaskName"]
                        if !occursin(task_name, json_file)
                            push!(warnings, "⚠️ $json_file: TaskName '$task_name' no coincide con el nombre del archivo")
                        end
                    end
                    
                catch e
                    push!(errors, "❌ Error parseando $json_file: $e")
                end
            end
        end
    end
    
    # Verificación opcional de archivos de tareas en raíz (convención de este proyecto).
    task_files = ["task-eyesopen_eeg.json", "task-eyesclosed_eeg.json"]
    for task_file in task_files
        if isfile(joinpath(bids_dir, task_file))
            println("✅ $task_file encontrado")
        else
            push!(warnings, "⚠️ Archivo de tarea faltante: $task_file")
        end
    end
    
    return errors, warnings
end

function main()
    # Directorio esperado por defecto cuando se ejecuta como script.
    bids_dir = "bids"
    
    if !isdir(bids_dir)
        println("❌ Directorio BIDS no encontrado: $bids_dir")
        return
    end
    
    println("🚀 Iniciando validación BIDS para: $bids_dir")
    println("="^60)
    
    # Ejecuta la validación y centraliza resultados para el reporte final.
    errors, warnings = validate_bids_structure(bids_dir)
    
    println("\n" * "="^60)
    println("📊 RESUMEN DE VALIDACIÓN")
    println("="^60)
    
    if isempty(errors)
        println("✅ ¡No se encontraron errores críticos!")
    else
        println("❌ ERRORES CRÍTICOS ENCONTRADOS:")
        for error in errors
            println("   $error")
        end
    end
    
    if !isempty(warnings)
        println("\n⚠️ ADVERTENCIAS:")
        for warning in warnings
            println("   $warning")
        end
    end
    
    println("\n📈 ESTADÍSTICAS:")
    println("   • Errores críticos: $(length(errors))")
    println("   • Advertencias: $(length(warnings))")
    
    if isempty(errors)
        println("\n🎉 ¡El dataset BIDS es válido!")
    else
        println("\n🔧 Se requieren correcciones antes de usar el dataset.")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Evita ejecutar main() cuando el archivo se importa como módulo.
    main()
end
