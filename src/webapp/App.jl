# NeuroMIND/src/webapp/App.jl
# Dashboard web — Pipeline Navigator + Phase 3 QC Initial
#
# Uso:  launch_webapp(cfg; port=8080)

using Genie
using Genie.Router
using Genie.Renderer.Html: html
using Genie.Renderer.Json: json
using Genie.Requests: getpayload
import Genie.Server
using Base64, Dates

"""
    launch_webapp(cfg; port=8080, open_browser=true)

Lanza el dashboard web de NeuroMIND en http://localhost:port.
"""
function launch_webapp(cfg::PipelineConfig;
                       port::Int      = 8080,
                       open_browser::Bool = true)

    res_root     = results_dir(cfg)
    project_root = cfg.root
    bids_root    = joinpath(res_root, "subjects")
    html_path    = joinpath(project_root, "web", "views", "dashboard.html")

    # ── Helpers internos ──────────────────────────────────────

    function _normalize_cond(cond::String)::String
        d = Dict(
            "EC"         => "eyesclosed",
            "EO"         => "eyesopen",
            "ec"         => "eyesclosed",
            "eo"         => "eyesopen",
            "eyesclosed" => "eyesclosed",
            "eyesopen"   => "eyesopen",
        )
        get(d, cond, lowercase(cond))
    end

    function _serve_csv(path::String)
        if isfile(path)
            try
                df   = CSV.read(path, DataFrame)
                rows = [Dict(zip(names(df),
                               [v isa Missing ? nothing : (v isa Bool ? v : string(v))
                                for v in collect(r)]))
                        for r in eachrow(df)]
                json(Dict("ok" => true, "columns" => names(df), "rows" => rows))
            catch e
                json(Dict("ok" => false, "rows" => [], "columns" => [],
                          "error" => string(e)))
            end
        else
            json(Dict("ok" => false, "rows" => [], "columns" => [],
                      "missing" => basename(path)))
        end
    end

    # ─── Helpers internos (Phase 4) ──────────────────────────

    function _apply_filter_chain(sig::Vector{Float64}, fs::Float64)::Vector{Float64}
        f   = cfg.filtering
        nyq = fs / 2.0
        ord = Int(get(f, "filter_order", 4))
        hp  = Float64(get(f, "highpass_hz",    0.5))
        lp  = Float64(get(f, "lowpass_hz",    48.0))
        nz  = Float64(get(f, "notch_hz",      50.0))
        nbw = Float64(get(f, "notch_bw_hz",    2.0))
        lo  = Float64(get(f, "bandreject_lo", 0.0))
        hi  = Float64(get(f, "bandreject_hi", 0.0))
        out = copy(sig)
        out = filtfilt(digitalfilter(Highpass(hp/nyq), Butterworth(ord)), out)
        out = filtfilt(digitalfilter(Lowpass(lp/nyq),  Butterworth(ord)), out)
        nz > 0 && (out = filtfilt(digitalfilter(
            Bandstop((nz - nbw/2)/nyq, (nz + nbw/2)/nyq), Butterworth(2)), out))
        lo > 0 && hi > lo && (out = filtfilt(digitalfilter(
            Bandstop(lo/nyq, hi/nyq), Butterworth(ord)), out))
        return out
    end

    function _welch_psd(sig::Vector{Float64}, fs::Float64, nfft::Int)
        nfft > length(sig) && return Float64[], Float64[]
        win     = Float64[0.5 - 0.5*cos(2π*(i-1)/(nfft-1)) for i in 1:nfft]
        win_pow = sum(win .^ 2)
        nhalf   = div(nfft, 2) + 1
        ps      = zeros(nhalf)
        n_seg   = div(length(sig), nfft)
        n_seg == 0 && return Float64[], Float64[]
        for i in 1:n_seg
            seg = sig[(i-1)*nfft+1 : i*nfft] .* win
            S   = abs.(rfft(seg)) .^ 2
            ps .+= S
        end
        ps ./= (n_seg * win_pow * fs)
        ps[2:end-1] .*= 2
        freqs = Float64[(k-1)*fs/nfft for k in 1:nhalf]
        return freqs, ps
    end

    function _read_channel_from_tsv(tsv_path::String, ch_name::String, n_pts::Int)
        sig = Float64[]
        open(tsv_path) do f
            readline(f)
            for line in eachline(f)
                isempty(strip(line)) && continue
                parts = split(line, '\t')
                length(parts) < 2 && continue
                string(strip(parts[1])) == ch_name || continue
                lim = min(n_pts, length(parts) - 1)
                sizehint!(sig, lim)
                for i in 1:lim
                    push!(sig, parse(Float64, parts[i+1]))
                end
                break
            end
        end
        return sig
    end

    # ─── API: configuración de filtros ────────────────────────
    route("/api/filter_config") do
        flt = cfg.filtering
        hp  = Float64(get(flt, "highpass_hz",    0.5))
        lp  = Float64(get(flt, "lowpass_hz",    48.0))
        nz  = Float64(get(flt, "notch_hz",      50.0))
        nbw = Float64(get(flt, "notch_bw_hz",    2.0))
        lo  = Float64(get(flt, "bandreject_lo", 100.0))
        hi  = Float64(get(flt, "bandreject_hi", 120.0))
        ord = Int(get(flt, "filter_order", 4))
        filters = [
            Dict("name"=>"High-pass",    "type"=>"Butterworth",
                 "freq"=>"$(hp) Hz",                    "order"=>ord, "applied"=>true),
            Dict("name"=>"Low-pass",     "type"=>"Butterworth",
                 "freq"=>"$(lp) Hz",                    "order"=>ord, "applied"=>true),
            Dict("name"=>"Notch 50 Hz",  "type"=>"IIR Notch",
                 "freq"=>"$(nz-nbw/2) – $(nz+nbw/2) Hz","order"=>2,   "applied"=>nz>0),
            Dict("name"=>"Notch 100 Hz", "type"=>"IIR Notch",
                 "freq"=>"$(lo) – $(hi) Hz",            "order"=>ord, "applied"=>lo>0&&hi>lo),
        ]
        json(Dict("ok"=>true, "filters"=>filters,
                  "highpass_hz"=>hp, "lowpass_hz"=>lp,
                  "notch_hz"=>nz,   "notch_bw_hz"=>nbw,
                  "bandreject_lo"=>lo, "bandreject_hi"=>hi,
                  "filter_order"=>ord))
    end

    # ─── API: señal cruda + filtrada de un canal ──────────────
    route("/api/bids/channel_signal") do
        subj    = string(get(getpayload(), :subj,    "M05"))
        sess    = string(get(getpayload(), :sess,    "T2"))
        cond    = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        ch_name = string(get(getpayload(), :channel, "Cz"))
        n_secs  = min(30, max(1, parse(Int, string(get(getpayload(), :secs, "10")))))

        tsv_name = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01_eeg_data.tsv"
        tsv_path = joinpath(project_root, "data", "BIDS", "raw", tsv_name)
        !isfile(tsv_path) && return json(Dict("ok"=>false,
            "error"=>"TSV no encontrado: $tsv_name"))

        try
            fs    = Float64(get(cfg.recording, "sampling_rate", 500.0))
            n_pts = Int(round(fs * n_secs))
            step  = max(1, div(n_pts, 900))

            raw = _read_channel_from_tsv(tsv_path, ch_name, n_pts)
            isempty(raw) && return json(Dict("ok"=>false,
                "error"=>"Canal $ch_name no encontrado en el TSV"))

            filt     = _apply_filter_chain(raw, fs)
            idx      = 1:step:length(raw)
            time_vec = Float64[(i-1)*step/fs for i in 1:length(idx)]

            json(Dict("ok"=>true, "channel"=>ch_name, "n_secs"=>n_secs,
                      "time_s"   => time_vec,
                      "raw"      => Float64[raw[i]  for i in idx],
                      "filtered" => Float64[filt[i] for i in idx],
                      "fs"       => fs))
        catch e
            json(Dict("ok"=>false, "error"=>string(e)))
        end
    end

    # ─── API: PSD cruda + filtrada de un canal ────────────────
    route("/api/bids/channel_psd") do
        subj    = string(get(getpayload(), :subj,    "M05"))
        sess    = string(get(getpayload(), :sess,    "T2"))
        cond    = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        ch_name = string(get(getpayload(), :channel, "Cz"))

        try
            fs   = Float64(get(cfg.recording, "sampling_rate", 500.0))
            nfft = Int(get(cfg.spectral, "nfft", 1024))

            # PSD filtrada desde CSV precomputado
            psd_path   = joinpath(bids_root, "sub-$subj", "ses-$sess", cond, "psd_by_channel.csv")
            filt_freqs = Float64[]; filt_power = Float64[]
            if isfile(psd_path)
                df   = CSV.read(psd_path, DataFrame)
                mask = df.channel .== ch_name
                filt_freqs = Float64.(df.freq_hz[mask])
                filt_power = Float64.(df.power_uv2[mask])
            end

            # PSD cruda calculada al vuelo
            tsv_name  = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01_eeg_data.tsv"
            tsv_path  = joinpath(project_root, "data", "BIDS", "raw", tsv_name)
            raw_freqs = Float64[]; raw_power = Float64[]
            band_impact = Dict{String,Any}()

            if isfile(tsv_path)
                raw = _read_channel_from_tsv(tsv_path, ch_name, 60_000)
                if !isempty(raw)
                    raw_freqs, raw_ps = _welch_psd(raw, fs, nfft)
                    raw_power         = raw_ps

                    # PSD de la señal filtrada para calcular impacto
                    filt_full = _apply_filter_chain(raw, fs)
                    _, flt_ps = _welch_psd(filt_full, fs, nfft)

                    bpw(freqs, ps, f1, f2) =
                        sum(ps[i] for (i, f) in enumerate(freqs) if f1 <= f <= f2; init=0.0)

                    band_impact = Dict(
                        "below_hp"     => round(100*(1 - bpw(raw_freqs,flt_ps,0,0.5)   /
                                                       max(1e-12, bpw(raw_freqs,raw_ps,0,0.5))),   digits=1),
                        "around_notch" => round(100*(1 - bpw(raw_freqs,flt_ps,48,52)   /
                                                       max(1e-12, bpw(raw_freqs,raw_ps,48,52))),   digits=1),
                        "above_lp"     => round(100*(1 - bpw(raw_freqs,flt_ps,48,fs/2) /
                                                       max(1e-12, bpw(raw_freqs,raw_ps,48,fs/2))), digits=1),
                        "signal_band"  => round(100* bpw(raw_freqs,flt_ps,1,40)        /
                                                       max(1e-12, bpw(raw_freqs,raw_ps,1,40)),     digits=1),
                    )
                end
            end

            json(Dict("ok"=>true, "channel"=>ch_name,
                      "raw_freqs"   => raw_freqs,  "raw_power"   => raw_power,
                      "filt_freqs"  => filt_freqs, "filt_power"  => filt_power,
                      "band_impact" => band_impact))
        catch e
            json(Dict("ok"=>false, "error"=>string(e)))
        end
    end

    # ─── API: potencia por banda desde CSV ────────────────────
    route("/api/bids/band_power") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        path = joinpath(bids_root, "sub-$subj", "ses-$sess", cond, "band_power_summary.csv")
        _serve_csv(path)
    end

    # ─── Dashboard HTML ───────────────────────────────────────
    route("/") do
        if isfile(html_path)
            html(read(html_path, String); noparse=true)
        else
            html("""
            <html><body style="font-family:sans-serif;padding:40px;background:#f0f4f8">
            <h2>⚠ dashboard.html no encontrado</h2>
            <p>Ruta esperada: <code>$(html_path)</code></p>
            <p>Asegúrate de que el archivo <code>web/views/dashboard.html</code> existe en el proyecto.</p>
            </body></html>"""; noparse=true)
        end
    end

    # ─── Config API ───────────────────────────────────────────
    route("/api/config") do
        bands_list = sort(collect(keys(cfg.bands)))
        json(Dict(
            "bands"        => bands_list,
            "qc_threshold" => 3.0,
            "version"      => "0.2.0",
        ))
    end

    # ─── BIDS API: sujetos ────────────────────────────────────
    route("/api/bids/subjects") do
        subjects = String[]
        if isdir(bids_root)
            for d in sort(readdir(bids_root))
                if startswith(d, "sub-") && isdir(joinpath(bids_root, d))
                    push!(subjects, d[5:end])
                end
            end
        end
        json(Dict("subjects" => subjects))
    end

    # ─── BIDS API: sesiones ───────────────────────────────────
    route("/api/bids/sessions") do
        subj     = string(get(getpayload(), :subj, ""))
        subj_dir = joinpath(bids_root, "sub-$(subj)")
        sessions = String[]
        if isdir(subj_dir)
            for d in sort(readdir(subj_dir))
                if startswith(d, "ses-") && isdir(joinpath(subj_dir, d))
                    push!(sessions, d[5:end])
                end
            end
        end
        json(Dict("sessions" => sessions))
    end

    # ─── BIDS API: condiciones ────────────────────────────────
    route("/api/bids/conditions") do
        subj     = string(get(getpayload(), :subj, ""))
        sess     = string(get(getpayload(), :sess, ""))
        sess_dir = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)")
        conds    = String[]
        if isdir(sess_dir)
            for d in sort(readdir(sess_dir))
                if isdir(joinpath(sess_dir, d))
                    push!(conds, d)
                end
            end
        end
        json(Dict("conditions" => conds))
    end

    # ─── BIDS API: estadísticas de canal ──────────────────────
    route("/api/bids/channel_stats") do
        subj = string(get(getpayload(), :subj, ""))
        sess = string(get(getpayload(), :sess, ""))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, "channel_statistics.csv")
        _serve_csv(path)
    end

    # ─── BIDS API: overview ───────────────────────────────────
    route("/api/bids/overview") do
        subj = string(get(getpayload(), :subj, ""))
        sess = string(get(getpayload(), :sess, ""))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, "overview.csv")
        _serve_csv(path)
    end

    # ─── BIDS API: log ────────────────────────────────────────
    route("/api/bids/log") do
        subj = string(get(getpayload(), :subj, ""))
        sess = string(get(getpayload(), :sess, ""))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, "pipeline_log.txt")
        if isfile(path)
            json(Dict("ok" => true, "content" => read(path, String)))
        else
            json(Dict("ok" => false, "content" => "Log no disponible para este sujeto/sesión."))
        end
    end

    # ─── BIDS API: señal cruda (TSV transpuesto) ─────────────
    # Formato TSV: filas=canales, columnas=muestras (T1, T2, ...)
    route("/api/bids/raw_signal") do
        subj   = string(get(getpayload(), :subj, "M05"))
        sess   = string(get(getpayload(), :sess, "T2"))
        cond   = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        n_secs = min(30, max(1, parse(Int, string(get(getpayload(), :secs, "10")))))

        task     = cond   # "eyesclosed" | "eyesopen"
        tsv_name = "sub-$(subj)_ses-$(sess)_task-$(task)_run-01_eeg_data.tsv"
        tsv_path = joinpath(project_root, "data", "BIDS", "raw", tsv_name)

        if !isfile(tsv_path)
            return json(Dict("ok" => false,
                             "error" => "Archivo no encontrado: $(tsv_name)"))
        end

        try
            fs    = Float64(get(cfg.recording, "sampling_rate", 500.0))
            n_pts = Int(round(fs * n_secs))   # muestras a leer
            step  = max(1, div(n_pts, 1000))  # decimación → ≤1000 puntos display

            ch_names = String[]
            matrix   = Vector{Float64}[]

            open(tsv_path) do f
                readline(f)  # saltar cabecera (Channel, T1, T2, ...)
                for line in eachline(f)
                    isempty(strip(line)) && continue
                    parts = split(line, '\t')
                    length(parts) < 2 && continue
                    ch = string(strip(parts[1]))
                    # Leer las primeras n_pts muestras y decimarlas
                    vals = Float64[]
                    limit = min(n_pts, length(parts) - 1)
                    for i in 1:step:limit
                        push!(vals, parse(Float64, parts[i + 1]))
                    end
                    push!(ch_names, ch)
                    push!(matrix, vals)
                end
            end

            n_disp   = isempty(matrix) ? 0 : length(matrix[1])
            time_vec = Float64[(i - 1) * step / fs for i in 1:n_disp]

            json(Dict(
                "ok"       => true,
                "channels" => ch_names,
                "time_s"   => time_vec,
                "matrix"   => matrix,
                "fs"       => fs,
                "n_secs"   => n_secs,
                "step"     => step,
            ))
        catch e
            json(Dict("ok" => false, "error" => string(e)))
        end
    end

    # ─── BIDS API: imagen (base64) ────────────────────────────
    route("/api/bids/image") do
        subj = string(get(getpayload(), :subj, ""))
        sess = string(get(getpayload(), :sess, ""))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        file = string(get(getpayload(), :file, ""))
        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, file)
        if isfile(path) && any(endswith(path, ext) for ext in [".png", ".jpg", ".svg"])
            json(Dict("ok" => true, "data" => base64encode(read(path))))
        else
            json(Dict("ok" => false, "error" => "No encontrado: $(file)"))
        end
    end

    # ─── BIDS API: estado del pipeline ────────────────────────
    route("/api/bids/pipeline_status") do
        subj = string(get(getpayload(), :subj, ""))
        sess = string(get(getpayload(), :sess, ""))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        base = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)

        fe(f) = isfile(joinpath(base, f))

        n_bad  = 0
        bad_ch = ""
        dur_s  = 0.0
        n_ch   = 0
        if fe("overview.csv")
            try
                df    = CSV.read(joinpath(base, "overview.csv"), DataFrame)
                n_bad = Int(df[1, :n_bad_ch])
                bad_ch = string(df[1, :bad_channels])
                dur_s  = Float64(df[1, :duration_s])
                n_ch   = Int(df[1, :n_channels])
            catch; end
        end

        # Inferir estado de cada fase a partir de archivos existentes
        statuses = Dict{String,String}(
            "0"  => isfile(joinpath(project_root, "config", "single_subject.toml")) ?
                    "completed" : "pending",
            "1"  => fe("overview.csv") ? "completed" : "pending",
            "2"  => fe("overview.csv") ? "completed" : "pending",
            "3"  => fe("channel_statistics.csv") ?
                    (n_bad > 0 ? "warning" : "completed") : "pending",
            "4"  => fe("filtered_signal_preview.png") ? "completed" : "pending",
            "5"  => "not_available",
            "6"  => fe("pipeline_log.txt") ? "completed" : "pending",
            "7"  => fe("pipeline_log.txt") ? "completed" : "pending",
            "8"  => fe("psd_by_channel.csv") ? "completed" : "pending",
            "9"  => fe("wpli_ALPHA.csv") ? "completed" : "pending",
            "10" => "not_available",
            "11" => fe("band_power_summary.csv") ? "completed" : "pending",
            "12" => "pending",
        )

        json(Dict(
            "ok"           => true,
            "statuses"     => statuses,
            "n_bad"        => n_bad,
            "bad_channels" => bad_ch,
            "duration_s"   => dur_s,
            "n_channels"   => n_ch,
        ))
    end

    # ─── API: BIDS & Metadata (Fase 1) ──────────────────────────
    route("/api/bids_metadata") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        raw_dir  = joinpath(project_root, "data", "BIDS", "raw")
        elec_dir = joinpath(project_root, "data", "BIDS", "electrodes")
        prefix   = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01"
        tsv_file  = "$(prefix)_eeg_data.tsv"
        meta_file = "$(prefix)_metadata.json"
        elec_file = "sub-$(subj)_ses-$(sess)_electrodes.tsv"

        tsv_path  = joinpath(raw_dir,  tsv_file)
        meta_path = joinpath(raw_dir,  meta_file)
        elec_path = joinpath(elec_dir, elec_file)

        # Metadata como raw string (se parsea en JS con JSON.parse)
        meta_raw = isfile(meta_path) ? read(meta_path, String) : ""

        # Electrodos TSV → array de dicts
        electrodes = Dict{String,Any}[]
        if isfile(elec_path)
            try
                df = CSV.read(elec_path, DataFrame; delim='\t')
                for r in eachrow(df)
                    push!(electrodes, Dict(
                        "name" => string(get(r, :name,  "")),
                        "x"    => r[:x] isa Missing ? nothing : Float64(r[:x]),
                        "y"    => r[:y] isa Missing ? nothing : Float64(r[:y]),
                        "z"    => r[:z] isa Missing ? nothing : Float64(r[:z]),
                        "type" => r[:type] isa Missing ? "EEG" : string(r[:type]),
                    ))
                end
            catch e
                @warn "Error reading electrodes TSV: $e"
            end
        end

        # Log del pipeline para timestamps
        res_base = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        log_path = joinpath(res_base, "pipeline_log.txt")
        log_content = isfile(log_path) ? read(log_path, String) : ""

        # Extraer inicio y duración del log
        log_start = ""; log_dur = ""
        if !isempty(log_content)
            m = match(r"NeuroMIND pipeline — (\S+)", log_content)
            m !== nothing && (log_start = string(m.captures[1]))
            m2 = match(r"completado en ([\d.]+) s", log_content)
            m2 !== nothing && (log_dur = string(m2.captures[1]) * " s")
        end

        # Tamaños de archivo
        sz(p) = isfile(p) ? round(filesize(p)/1024, digits=1) : 0.0

        json(Dict(
            "ok"        => true,
            "meta_json" => meta_raw,
            "electrodes"=> electrodes,
            "file_tree" => Dict(
                "subj"       => "sub-$(subj)",
                "sess"       => "ses-$(sess)",
                "tsv"        => Dict("name"=>tsv_file,  "exists"=>isfile(tsv_path),  "size_kb"=>sz(tsv_path)),
                "metadata"   => Dict("name"=>meta_file, "exists"=>isfile(meta_path), "size_kb"=>sz(meta_path)),
                "electrodes" => Dict("name"=>elec_file, "exists"=>isfile(elec_path), "size_kb"=>sz(elec_path)),
            ),
            "log_start" => log_start,
            "log_dur"   => log_dur,
            "has_log"   => isfile(log_path),
        ))
    end

    # ─── API: resumen del proyecto (Fase 0) ──────────────────────
    route("/api/project_overview") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        raw_dir  = joinpath(project_root, "data", "BIDS", "raw")
        elec_dir = joinpath(project_root, "data", "BIDS", "electrodes")
        prefix   = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01"
        tsv_file = "$(prefix)_eeg_data.tsv"
        meta_file= "$(prefix)_metadata.json"
        elec_file= "sub-$(subj)_ses-$(sess)_electrodes.tsv"

        # Metadata JSON como string crudo (parseado en JS con JSON.parse)
        meta_raw = ""
        meta_path = joinpath(raw_dir, meta_file)
        isfile(meta_path) && (meta_raw = read(meta_path, String))

        # Disponibilidad de resultados procesados
        res_base = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        fe(f)    = isfile(joinpath(res_base, f))

        # Sujetos y sesiones con resultados
        available_subj = String[]
        if isdir(bids_root)
            for d in sort(readdir(bids_root))
                startswith(d, "sub-") && isdir(joinpath(bids_root, d)) &&
                    push!(available_subj, d)
            end
        end
        n_sessions = 0
        for sd in available_subj
            sdir = joinpath(bids_root, sd)
            isdir(sdir) || continue
            for d in readdir(sdir)
                startswith(d, "ses-") && isdir(joinpath(sdir, d)) && (n_sessions += 1)
            end
        end

        # Sujetos procesados desde subjects_index.csv
        n_processed = 0
        idx_path = joinpath(res_root, "subjects_index.csv")
        if isfile(idx_path)
            try
                df = CSV.read(idx_path, DataFrame)
                n_processed = nrow(df)
            catch; end
        end

        json(Dict(
            "ok"               => true,
            "meta_json"        => meta_raw,
            "file_tree"        => Dict(
                "tsv"        => Dict("name"=>tsv_file,  "exists"=>isfile(joinpath(raw_dir,  tsv_file))),
                "metadata"   => Dict("name"=>meta_file, "exists"=>isfile(joinpath(raw_dir,  meta_file))),
                "electrodes" => Dict("name"=>elec_file, "exists"=>isfile(joinpath(elec_dir, elec_file))),
            ),
            "results"          => Dict(
                "overview"      => fe("overview.csv"),
                "channel_stats" => fe("channel_statistics.csv"),
                "psd"           => fe("psd_by_channel.csv"),
                "wpli"          => fe("wpli_ALPHA.csv"),
                "band_power"    => fe("band_power_summary.csv"),
            ),
            "n_subjects_bids"  => length(available_subj),
            "n_sessions_total" => n_sessions,
            "n_processed"      => n_processed,
        ))
    end

    # ─── API: información señal cruda (Fase 2) ──────────────────
    route("/api/phase2_info") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        raw_dir   = joinpath(project_root, "data", "BIDS", "raw")
        prefix    = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01"
        tsv_file  = "$(prefix)_eeg_data.tsv"
        meta_file = "$(prefix)_metadata.json"
        tsv_path  = joinpath(raw_dir, tsv_file)
        meta_path = joinpath(raw_dir, meta_file)

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        ov_path   = joinpath(res_base, "overview.csv")
        log_path  = joinpath(res_base, "pipeline_log.txt")

        # ── Signal info desde overview.csv o metadata ──────────
        fs = Float64(get(cfg.recording, "sampling_rate", 500.0))
        n_samples = 0; n_channels = 0; duration_s = 0.0
        n_bad = 0; bad_str = ""

        if isfile(ov_path)
            try
                df = CSV.read(ov_path, DataFrame)
                if nrow(df) > 0
                    r = df[1, :]
                    fs         = Float64(r.sampling_hz)
                    n_channels = Int(r.n_channels)
                    n_samples  = Int(r.n_samples)
                    duration_s = Float64(r.duration_s)
                    n_bad      = Int(r.n_bad_ch)
                    bad_str    = string(get(r, :bad_channels, ""))
                end
            catch; end
        elseif isfile(meta_path)
            try
                txt = read(meta_path, String)
                for (pat, sym) in [
                    (r"\"fs\"\s*:\s*([\d.]+)",        :fs),
                    (r"\"n_samples\"\s*:\s*(\d+)",     :n_samp),
                    (r"\"n_channels\"\s*:\s*(\d+)",    :n_ch),
                    (r"\"duration_s\"\s*:\s*([\d.]+)", :dur),
                ]
                    m = match(pat, txt); m === nothing && continue
                    val = m.captures[1]
                    sym === :fs     && (fs         = parse(Float64, val))
                    sym === :n_samp && (n_samples   = parse(Int, val))
                    sym === :n_ch   && (n_channels  = parse(Int, val))
                    sym === :dur    && (duration_s  = parse(Float64, val))
                end
            catch; end
        end

        # ── Estadísticas de amplitud: muestrea ~10 s crudos ────
        g_min = Inf; g_max = -Inf
        g_sum = 0.0; g_sum2 = 0.0; g_count = 0
        sat_ch = 0
        if isfile(tsv_path)
            max_pts = Int(round(fs * 10.0))
            try
                open(tsv_path) do f
                    readline(f)
                    for line in eachline(f)
                        isempty(strip(line)) && continue
                        parts = split(line, '\t')
                        length(parts) < 2 && continue
                        lim = min(max_pts, length(parts) - 1)
                        ch_lo = Inf; ch_hi = -Inf
                        for i in 1:lim
                            v = parse(Float64, strip(parts[i + 1]))
                            v < ch_lo && (ch_lo = v)
                            v > ch_hi && (ch_hi = v)
                            g_sum  += v
                            g_sum2 += v * v
                            g_count += 1
                        end
                        ch_lo < g_min && (g_min = ch_lo)
                        ch_hi > g_max && (g_max = ch_hi)
                        (ch_hi - ch_lo) > 800 && (sat_ch += 1)
                    end
                end
            catch; end
        end

        g_min  = isinf(g_min) ? 0.0 : g_min
        g_max  = isinf(g_max) ? 0.0 : g_max
        g_mean = g_count > 0 ? g_sum  / g_count : 0.0
        g_rms  = g_count > 0 ? sqrt(g_sum2 / g_count) : 0.0
        g_p2p  = g_max - g_min

        # ── Duración formateada ──
        dur_fmt = ""
        if duration_s > 0
            mm = div(Int(floor(duration_s)), 60)
            ss = Int(round(duration_s)) % 60
            dur_fmt = "$(mm):$(lpad(ss, 2, '0')) min ($(Int(round(duration_s))) s)"
        end

        # ── Timing de fase ──
        log_start = ""; log_dur_s = ""
        log_ok = isfile(log_path)
        if log_ok
            txt = read(log_path, String)
            m  = match(r"NeuroMIND pipeline — (\S+)", txt)
            m  !== nothing && (log_start = string(m.captures[1]))
            m2 = match(r"completado en ([\d.]+) s", txt)
            m2 !== nothing && (log_dur_s = string(m2.captures[1]))
        end

        tsv_exists  = isfile(tsv_path)
        tsv_size_mb = tsv_exists ? round(filesize(tsv_path) / 1048576, digits=1) : 0.0

        json(Dict(
            "ok" => true,
            "signal_info" => Dict(
                "duration_s"   => round(duration_s, digits=2),
                "duration_fmt" => dur_fmt,
                "fs_hz"        => fs,
                "n_samples"    => n_samples,
                "n_channels"   => n_channels,
                "reference"    => "Cz",
                "montage"      => "custom",
                "unit"         => "µV",
            ),
            "amplitude" => Dict(
                "global_min"   => round(g_min,  digits=1),
                "global_max"   => round(g_max,  digits=1),
                "peak_to_peak" => round(g_p2p,  digits=1),
                "rms_global"   => round(g_rms,  digits=1),
                "mean_global"  => round(g_mean, digits=1),
            ),
            "status" => Dict(
                "loaded"           => tsv_exists,
                "no_saturation"    => sat_ch == 0,
                "no_interruptions" => true,
                "quality"          => n_bad == 0 ? "good" : "warning",
                "n_bad_channels"   => n_bad,
                "bad_channels"     => isempty(bad_str) ?
                                      String[] :
                                      [strip(b) for b in split(bad_str, ",")],
            ),
            "file_info" => Dict(
                "filename" => tsv_file,
                "format"   => "TSV (BIDS)",
                "size_mb"  => tsv_size_mb,
                "exists"   => tsv_exists,
            ),
            "phase_timing" => Dict(
                "start"    => log_start,
                "duration" => log_dur_s,
                "has_log"  => log_ok,
            ),
        ))
    end

    # ─── API: Fase 5 — ICA ──────────────────────────────────────
    route("/api/phase5_ica_info") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        comp_path = joinpath(res_base, "ica_components.csv")
        summ_path = joinpath(res_base, "ica_summary.json")
        log_path  = joinpath(res_base, "pipeline_log.txt")

        n_comp_cfg = Int(get(cfg.ica, "n_components", 30))
        method_cfg = string(get(cfg.ica, "method", "fastica"))

        ica_config = Dict(
            "n_components" => n_comp_cfg,
            "method"       => uppercase(method_cfg),
            "algorithm"    => "PCA whitening + FastICA",
            "library"      => "LinearAlgebra (Julia puro)",
            "max_iter"     => 500,
            "tol"          => 1e-5,
            "seed"         => 42,
        )

        if !isfile(comp_path)
            return json(Dict(
                "ok"          => true,
                "ica_run"     => false,
                "ica_config"  => ica_config,
                "components"  => Dict{String,Any}[],
                "summary"     => Dict(
                    "n_components"      => 0,
                    "n_rejected"        => 0,
                    "n_accepted"        => 0,
                    "variance_retained" => 0.0,
                    "run_timestamp"     => "",
                    "run_duration_s"    => 0.0,
                ),
                "figures"      => String[],
                "phase_timing" => Dict("start"=>"", "duration"=>0.0, "has_log"=>false),
            ))
        end

        components = Dict{String,Any}[]
        try
            df = CSV.read(comp_path, DataFrame)
            for row in eachrow(df)
                push!(components, Dict{String,Any}(
                    "index"        => Int(get(row, :component, get(row, :index, 0))),
                    "label"        => string(get(row, :label, "IC")),
                    "type"         => string(get(row, :artifact_type, get(row, :type, "unknown"))),
                    "variance_pct" => round(Float64(get(row, :variance_pct,
                                                        get(row, :variance, 0.0))), digits=2),
                    "rejected"     => Bool(get(row, :rejected, false)),
                ))
            end
        catch e
            @warn "ICA components CSV parse error: $e"
        end

        n_rej   = count(c -> Bool(get(c, "rejected", false)), components)
        n_acc   = length(components) - n_rej
        var_rej = sum(c -> Bool(get(c, "rejected", false)) ?
                             Float64(get(c, "variance_pct", 0.0)) : 0.0,
                     components; init=0.0)
        var_ret = round(100.0 - var_rej, digits=1)
        run_ts  = ""; run_dur = 0.0
        if isfile(summ_path)
            try
                txt   = read(summ_path, String)
                m_ts  = match(r"\"timestamp\"\s*:\s*\"([^\"]+)\"", txt)
                m_dur = match(r"\"duration_s\"\s*:\s*([0-9.eE+\-]+)", txt)
                if m_ts  !== nothing; run_ts  = String(m_ts.captures[1]); end
                if m_dur !== nothing; run_dur = parse(Float64, m_dur.captures[1]); end
            catch; end
        end

        figs_dir = joinpath(res_base, "figures")
        ica_figs = String[]
        if isdir(figs_dir)
            ica_figs = filter(
                f -> occursin(r"ica|ICA|component|IC"i, f) &&
                     (endswith(f, ".png") || endswith(f, ".svg")),
                readdir(figs_dir)
            )
        end

        log_start = ""; log_dur_s = 0.0; log_ok = false
        if isfile(log_path)
            try
                txt = read(log_path, String)
                m   = match(r"ICA.*?(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})", txt)
                if m !== nothing; log_start = m.captures[1]; log_ok = true; end
                m2  = match(r"ICA.*?(\d+\.?\d*)\s*s", txt)
                if m2 !== nothing; log_dur_s = parse(Float64, m2.captures[1]); end
            catch; end
        end

        json(Dict(
            "ok"      => true,
            "ica_run" => true,
            "ica_config"  => ica_config,
            "components"  => components,
            "summary" => Dict(
                "n_components"      => length(components),
                "n_rejected"        => n_rej,
                "n_accepted"        => n_acc,
                "variance_retained" => var_ret,
                "run_timestamp"     => run_ts,
                "run_duration_s"    => run_dur,
            ),
            "figures"      => ica_figs,
            "phase_timing" => Dict("start"=>log_start, "duration"=>log_dur_s, "has_log"=>log_ok),
        ))
    end

    # ─── API: Señal de componente ICA ───────────────────────────
    route("/api/ica_activation") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        comp = Int(parse(Float64, string(get(getpayload(), :comp, "1"))))

        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, "ica_activations.csv")
        if !isfile(path)
            return json(Dict("ok"=>false, "msg"=>"Sin activaciones guardadas"))
        end
        df  = CSV.read(path, DataFrame)
        col = Symbol("IC$(comp)")
        if !hasproperty(df, col)
            return json(Dict("ok"=>false, "msg"=>"Componente $(comp) no encontrado"))
        end
        t_vec  = Float64.(df.t_s)
        values = Float64.(df[!, col])
        json(Dict("ok"=>true, "t"=>t_vec, "values"=>values,
                  "fs"=>Float64(get(cfg.recording,"fs",500.0)),
                  "comp"=>comp))
    end

    route("/api/ica_signal") do
        subj  = string(get(getpayload(), :subj, "M05"))
        sess  = string(get(getpayload(), :sess, "T2"))
        cond  = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        stage = string(get(getpayload(), :stage, "before"))   # "before" | "after"
        n_ch  = Int(parse(Float64, string(get(getpayload(), :channels, "8"))))

        path  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond,
                         "ica_signal_$(stage).csv")
        if !isfile(path)
            return json(Dict("ok"=>false, "msg"=>"Sin señal $(stage) guardada"))
        end
        df = CSV.read(path, DataFrame)
        t_vec = Float64.(df.t_s)
        cols  = [n for n in names(df) if n != "t_s"][1:min(n_ch, ncol(df)-1)]
        channels = [Dict("name"=>c, "values"=>Float64.(df[!, c])) for c in cols]
        json(Dict("ok"=>true, "t"=>t_vec, "channels"=>channels))
    end

    # ─── Legacy API (mantener compatibilidad) ─────────────────

    route("/api/subjects") do
        dirs = isdir(res_root) ?
               filter(d -> isdir(joinpath(res_root, d)) &&
                            d ∉ ["group", "subjects"],
                      readdir(res_root)) : String[]
        json(Dict("subjects" => dirs))
    end

    route("/api/sessions/:subj") do
        subj_dir = joinpath(res_root, params(:subj))
        sess = isdir(subj_dir) ?
               filter(d -> isdir(joinpath(subj_dir, d)), readdir(subj_dir)) :
               String[]
        json(Dict("sessions" => sess))
    end

    route("/api/figures/:subj/:sess/:cond") do
        fig_dir = joinpath(res_root, params(:subj), params(:sess), "figures")
        cond    = params(:cond)
        figs    = isdir(fig_dir) ?
                  filter(f -> endswith(f, ".png") &&
                              (cond == "ALL" || isempty(cond) || occursin(cond, f)),
                         readdir(fig_dir)) : String[]
        json(Dict("figures" => figs))
    end

    route("/api/image") do
        rel      = string(get(getpayload(), :path, ""))
        abs_path = joinpath(res_root, rel)
        if isfile(abs_path) && endswith(abs_path, ".png")
            json(Dict("data" => base64encode(read(abs_path)), "ok" => true))
        else
            json(Dict("ok" => false, "error" => "No encontrado: $rel"))
        end
    end

    route("/api/table") do
        rel      = string(get(getpayload(), :path, ""))
        abs_path = joinpath(res_root, rel)
        if isfile(abs_path) && endswith(abs_path, ".csv")
            df   = CSV.read(abs_path, DataFrame)
            rows = [Dict(zip(names(df), collect(r))) for r in eachrow(df)]
            json(Dict("ok" => true, "columns" => names(df), "rows" => rows))
        else
            json(Dict("ok" => false, "rows" => [], "columns" => []))
        end
    end

    route("/api/export/:subj/:sess") do
        subj  = params(:subj)
        sess  = params(:sess)
        base  = joinpath(res_root, subj, sess)
        files = Dict{String,Vector{String}}(
            "tables"  => String[],
            "figures" => String[],
            "logs"    => String[],
        )
        for (key, subdir) in [("tables","tables"), ("figures","figures"), ("logs","logs")]
            d = joinpath(base, subdir)
            isdir(d) && (files[key] = readdir(d))
        end
        json(files)
    end

    route("/api/report/:subj/:sess/:cond") do
        subj     = params(:subj)
        sess     = params(:sess)
        cond     = params(:cond)
        rep_path = joinpath(res_root, subj, sess, "reports", "report_$(cond).html")
        if isfile(rep_path)
            json(Dict("ok" => true, "html" => read(rep_path, String)))
        else
            json(Dict("ok" => false, "html" => ""))
        end
    end

    # ─── Iniciar servidor ─────────────────────────────────────
    @info "NeuroMIND Dashboard → http://localhost:$(port)"
    open_browser && _try_open_browser("http://localhost:$(port)")
    Genie.Server.up(port; async=true)
    return nothing
end

function _try_open_browser(url::String)
    try
        Sys.isapple() && run(`open $url`, wait=false)
        Sys.islinux() && run(`xdg-open $url`, wait=false)
    catch; end
end
