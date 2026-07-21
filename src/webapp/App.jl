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
        f       = cfg.filtering
        nyq     = fs / 2.0
        profile = get(f, "profile", "default")
        ord     = Int(get(f, "filter_order", 4))
        hp      = Float64(get(f, "highpass_hz",    0.5))
        lp      = Float64(get(f, "lowpass_hz",   150.0))
        nz      = Float64(get(f, "notch_hz",      50.0))
        nbw     = Float64(get(f, "notch_bw_hz",    1.0))
        lo      = Float64(get(f, "bandreject_lo", 99.5))
        hi      = Float64(get(f, "bandreject_hi",100.5))
        out     = copy(sig)

        if profile == "eeg_julia"
            # EEG_Julia: Notch(filt) → Bandreject(filt) → HP(filtfilt) → LP(filtfilt)
            nz > 0 && (out = filt(digitalfilter(
                Bandstop((nz-nbw/2)/nyq, (nz+nbw/2)/nyq), Butterworth(ord)), out))
            lo > 0 && hi > lo && (out = filt(digitalfilter(
                Bandstop(lo/nyq, hi/nyq), Butterworth(ord)), out))
            out = filtfilt(digitalfilter(Highpass(hp/nyq), Butterworth(ord)), out)
            out = filtfilt(digitalfilter(Lowpass(lp/nyq),  Butterworth(ord)), out)
        else
            # Default: HP(filtfilt) → LP(filtfilt) → Notch(filtfilt) → Bandreject(filtfilt)
            out = filtfilt(digitalfilter(Highpass(hp/nyq), Butterworth(ord)), out)
            out = filtfilt(digitalfilter(Lowpass(lp/nyq),  Butterworth(ord)), out)
            nz > 0 && (out = filtfilt(digitalfilter(
                Bandstop((nz-nbw/2)/nyq, (nz+nbw/2)/nyq), Butterworth(ord)), out))
            lo > 0 && hi > lo && (out = filtfilt(digitalfilter(
                Bandstop(lo/nyq, hi/nyq), Butterworth(ord)), out))
        end
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
        flt     = cfg.filtering
        profile = get(flt, "profile", "default")
        hp      = Float64(get(flt, "highpass_hz",    0.5))
        lp      = Float64(get(flt, "lowpass_hz",   150.0))
        nz      = Float64(get(flt, "notch_hz",      50.0))
        nbw     = Float64(get(flt, "notch_bw_hz",    1.0))
        lo      = Float64(get(flt, "bandreject_lo",  99.5))
        hi      = Float64(get(flt, "bandreject_hi", 100.5))
        ord     = Int(get(flt, "filter_order", 4))

        # Cadena real según perfil (usa la misma lógica que filter_recording)
        chain = describe_filter_chain(cfg)
        filters = [Dict(
            "step"    => s.step,
            "name"    => s.name,
            "type"    => "Butterworth",
            "freq"    => s.freq,
            "order"   => s.order,
            "method"  => s.method,
            "applied" => true,
        ) for s in chain]

        json(Dict("ok"=>true,
                  "profile"        => profile,
                  "filters"        => filters,
                  "highpass_hz"    => hp,
                  "lowpass_hz"     => lp,
                  "notch_hz"       => nz,
                  "notch_bw_hz"    => nbw,
                  "bandreject_lo"  => lo,
                  "bandreject_hi"  => hi,
                  "filter_order"   => ord))
    end

    # ─── API: respuesta real de filtros (DSP.freqresp) ──────────
    # Devuelve magnitud en dB para cada filtro individual y la cadena compuesta.
    # Método: filt → |H(ω)|, filtfilt → |H(ω)|² (fase cero, orden efectivo doble).
    route("/api/filter_response") do
        flt_cfg = cfg.filtering
        profile = get(flt_cfg, "profile", "default")
        hp      = Float64(get(flt_cfg, "highpass_hz",    0.5))
        lp      = Float64(get(flt_cfg, "lowpass_hz",   150.0))
        nz      = Float64(get(flt_cfg, "notch_hz",      50.0))
        nbw     = Float64(get(flt_cfg, "notch_bw_hz",    1.0))
        lo      = Float64(get(flt_cfg, "bandreject_lo",  99.5))
        hi      = Float64(get(flt_cfg, "bandreject_hi", 100.5))
        ord     = Int(get(flt_cfg, "filter_order", 4))
        fs      = Float64(get(cfg.recording, "fs", 500.0))
        nyq     = fs / 2.0

        n_pts = 2048
        omega = Float64[π * k / (n_pts - 1) for k in 0:(n_pts - 1)]
        freqs = Float64[ω * nyq / π for ω in omega]

        f_notch = digitalfilter(Bandstop((nz - nbw/2)/nyq, (nz + nbw/2)/nyq), Butterworth(ord))
        f_br    = digitalfilter(Bandstop(lo/nyq, hi/nyq),                       Butterworth(ord))
        f_hp    = digitalfilter(Highpass(hp/nyq),                                Butterworth(ord))
        f_lp    = digitalfilter(Lowpass(lp/nyq),                                 Butterworth(ord))

        H_notch = freqresp(f_notch, omega)
        H_br    = freqresp(f_br,    omega)
        H_hp    = freqresp(f_hp,    omega)
        H_lp    = freqresp(f_lp,    omega)

        # filt (causal) → |H|; filtfilt (zero-phase) → |H|² (magnitude squared)
        m_notch = abs.(H_notch)
        m_br    = abs.(H_br)
        m_hp    = abs.(H_hp) .^ 2
        m_lp    = abs.(H_lp) .^ 2
        m_comp  = m_notch .* m_br .* m_hp .* m_lp

        to_db(m) = Float64[20.0 * log10(max(v, 1e-10)) for v in m]
        r4(v)    = round.(v, digits=4)

        json(Dict(
            "ok"      => true,
            "fs"      => fs,
            "freqs"   => r4(freqs),
            "profile" => profile,
            "individual" => Dict(
                "Notch"      => Dict(
                    "mag_db"    => r4(to_db(m_notch)),
                    "method"    => "filt",
                    "order_eff" => ord,
                    "freq_label"=> "$(nz - nbw/2)–$(nz + nbw/2) Hz",
                ),
                "Bandreject" => Dict(
                    "mag_db"    => r4(to_db(m_br)),
                    "method"    => "filt",
                    "order_eff" => ord,
                    "freq_label"=> "$(lo)–$(hi) Hz",
                ),
                "High-pass"  => Dict(
                    "mag_db"    => r4(to_db(m_hp)),
                    "method"    => "filtfilt",
                    "order_eff" => ord * 2,
                    "freq_label"=> "$(hp) Hz",
                ),
                "Low-pass"   => Dict(
                    "mag_db"    => r4(to_db(m_lp)),
                    "method"    => "filtfilt",
                    "order_eff" => ord * 2,
                    "freq_label"=> "$(lp) Hz",
                ),
            ),
            "composite" => Dict("mag_db" => r4(to_db(m_comp))),
            "markers"   => Dict(
                "notch_hz"      => nz,
                "notch_bw_hz"   => nbw,
                "bandreject_lo" => lo,
                "bandreject_hi" => hi,
                "highpass_hz"   => hp,
                "lowpass_hz"    => lp,
            ),
        ))
    end

    # ─── API: PSD por etapas de filtrado (estilo EEG_Julia) ─────
    # Aplica cada filtro secuencialmente y devuelve PSD (Welch-Hamming, nfft=N)
    # para cada etapa, replicando el método de EEG_Julia/src/Preprocessing/filtering.jl.
    route("/api/bids/channel_psd_stages") do
        subj    = string(get(getpayload(), :subj, "M05"))
        sess    = string(get(getpayload(), :sess, "T2"))
        cond    = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        ch_name = string(get(getpayload(), :channel, "Cz"))

        try
            fs  = Float64(get(cfg.recording, "fs", 500.0))
            nyq = fs / 2.0
            flt = cfg.filtering
            ord = Int(get(flt, "filter_order", 4))
            hp  = Float64(get(flt, "highpass_hz",    0.5))
            lp  = Float64(get(flt, "lowpass_hz",   150.0))
            nz  = Float64(get(flt, "notch_hz",      50.0))
            nbw = Float64(get(flt, "notch_bw_hz",    1.0))
            lo  = Float64(get(flt, "bandreject_lo",  99.5))
            hi  = Float64(get(flt, "bandreject_hi", 100.5))

            tsv_name = "sub-$(subj)_ses-$(sess)_task-$(cond)_run-01_eeg_data.tsv"
            tsv_path = joinpath(project_root, "data", "BIDS", "raw", tsv_name)
            !isfile(tsv_path) && return json(Dict("ok"=>false,
                "error"=>"TSV no encontrado: $tsv_name"))

            raw = _read_channel_from_tsv(tsv_path, ch_name, 60_000)
            isempty(raw) && return json(Dict("ok"=>false,
                "error"=>"Canal $ch_name no encontrado en $tsv_name"))

            # PSD estilo EEG_Julia: welch_pgram con ventana Hamming, nfft = N
            n_sig = length(raw)
            function _psd_stages(sig)
                p = welch_pgram(sig; fs=fs, window=hamming, nfft=n_sig)
                return Float64.(DSP.freq(p)), Float64.(DSP.power(p))
            end

            # Cadena "eeg_julia": Notch(filt) → BR(filt) → HP(filtfilt) → LP(filtfilt)
            s0 = raw
            s1 = filt(digitalfilter(Bandstop((nz-nbw/2)/nyq,(nz+nbw/2)/nyq),Butterworth(ord)), s0)
            s2 = filt(digitalfilter(Bandstop(lo/nyq, hi/nyq), Butterworth(ord)), s1)
            s3 = filtfilt(digitalfilter(Highpass(hp/nyq), Butterworth(ord)), s2)
            s4 = filtfilt(digitalfilter(Lowpass(lp/nyq),  Butterworth(ord)), s3)

            freqs, p0 = _psd_stages(s0)
            _,     p1 = _psd_stages(s1)
            _,     p2 = _psd_stages(s2)
            _,     p3 = _psd_stages(s3)
            _,     p4 = _psd_stages(s4)

            r5(v) = round.(v, digits=6)
            json(Dict(
                "ok"      => true,
                "channel" => ch_name,
                "method"  => "welch_pgram_hamming_nfft_N",
                "n_pts"   => n_sig,
                "freqs"   => r5(freqs),
                "stages"  => Dict(
                    "original"   => r5(p0),
                    "notch"      => r5(p1),
                    "bandreject" => r5(p2),
                    "highpass"   => r5(p3),
                    "lowpass"    => r5(p4),
                ),
                "stage_labels" => [
                    "Original",
                    "Notch 50 Hz",
                    "Notch + BR 100 Hz",
                    "Notch + BR + HP 0.5 Hz",
                    "Final (+ LP 150 Hz)",
                ],
            ))
        catch e
            json(Dict("ok"=>false, "error"=>string(e)))
        end
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
            "0"  => isfile(joinpath(project_root, "config", "pipeline.toml")) ?
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
            "10" => fe("surrogate_summary.json") ? "completed" :
                    fe("significant_connections.csv") ? "completed" : "pending",
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

        profile_cfg  = string(get(cfg.ica, "profile", "default"))
        n_comp_cfg   = Int(get(cfg.ica, "n_components", 30))
        method_cfg   = string(get(cfg.ica, "method", "fastica"))
        is_eeg_julia = profile_cfg == "eeg_julia"

        ica_config = Dict(
            "profile"      => profile_cfg,
            "n_components" => is_eeg_julia ? "n_channels" : n_comp_cfg,
            "method"       => uppercase(method_cfg),
            "algorithm"    => "PCA whitening + FastICA simétrico",
            "library"      => "LinearAlgebra (Julia puro)",
            "max_iter"     => is_eeg_julia ? 512 : Int(get(cfg.ica, "max_iter", 500)),
            "tol"          => is_eeg_julia ? 1e-7 : Float64(get(cfg.ica, "tol", 1e-5)),
            "seed"         => is_eeg_julia ? 1234 : Int(get(cfg.ica, "seed", 42)),
            "eeg_julia_compatible" => is_eeg_julia,
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

        # ── Componentes desde ica_components.csv ─────────────────
        components = Dict{String,Any}[]
        figs_dir   = joinpath(res_base, "figures")
        try
            df = CSV.read(comp_path, DataFrame)
            for row in eachrow(df)
                ic_idx = Int(get(row, :component, get(row, :index, 0)))
                # Buscar topomap correspondiente
                topo_name = "ica_topomap_$(lpad(ic_idx, 3, '0')).png"
                topo_url  = isfile(joinpath(figs_dir, topo_name)) ?
                            "/api/ica_topomap?subj=$(subj)&sess=$(sess)&cond=$(cond)&file=$(topo_name)" : ""
                push!(components, Dict{String,Any}(
                    "index"        => ic_idx,
                    "label"        => string(get(row, :label, "IC")),
                    "type"         => string(get(row, :artifact_type, get(row, :type, "unknown"))),
                    "variance_pct" => round(Float64(get(row, :variance_pct,
                                                        get(row, :variance, 0.0))), digits=2),
                    "rejected"     => Bool(get(row, :rejected, false)),
                    "topomap_url"  => topo_url,
                ))
            end
        catch e
            @warn "ICA components CSV parse error: $e"
        end

        # ── Métricas de resumen ───────────────────────────────
        n_rej   = count(c -> Bool(get(c, "rejected", false)), components)
        n_acc   = length(components) - n_rej
        var_rej = sum(c -> Bool(get(c, "rejected", false)) ?
                             Float64(get(c, "variance_pct", 0.0)) : 0.0,
                     components; init=0.0)
        var_ret     = round(100.0 - var_rej, digits=1)
        run_ts      = ""; run_dur = 0.0
        has_feat    = false; n_topo = 0; has_summ = false
        artifact_types_count = Dict{String,Int}()
        if isfile(summ_path)
            has_summ = true
            try
                txt   = read(summ_path, String)
                m_ts  = match(r"\"timestamp\"\s*:\s*\"([^\"]+)\"", txt)
                m_dur = match(r"\"duration_s\"\s*:\s*([0-9.eE+\-]+)", txt)
                m_feat= match(r"\"has_features\"\s*:\s*(true|false)", txt)
                m_top = match(r"\"n_topomaps\"\s*:\s*(\d+)", txt)
                if m_ts   !== nothing; run_ts  = String(m_ts.captures[1]); end
                if m_dur  !== nothing; run_dur = parse(Float64, m_dur.captures[1]); end
                if m_feat !== nothing; has_feat = m_feat.captures[1] == "true"; end
                if m_top  !== nothing; n_topo  = parse(Int, m_top.captures[1]); end
            catch; end
        end

        # Cuenta tipos de artefacto desde components
        for c in components
            t = string(get(c, "type", "unknown"))
            artifact_types_count[t] = get(artifact_types_count, t, 0) + 1
        end

        ica_figs = String[]
        if isdir(figs_dir)
            ica_figs = sort(filter(
                f -> occursin(r"ica|ICA|component|IC"i, f) &&
                     (endswith(f, ".png") || endswith(f, ".svg")),
                readdir(figs_dir)
            ))
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
                "n_components"       => length(components),
                "n_rejected"         => n_rej,
                "n_accepted"         => n_acc,
                "variance_retained"  => var_ret,
                "run_timestamp"      => run_ts,
                "run_duration_s"     => run_dur,
                "has_features"       => has_feat,
                "n_topomaps"         => n_topo,
                "artifact_types"     => artifact_types_count,
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

    # ─── API: Imagen de topomap ICA (base64 JSON) ────────────
    route("/api/ica_topomap") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        file = string(get(getpayload(), :file, ""))

        # Validación de seguridad: solo ica_topomap_NNN.png
        if !occursin(r"^ica_topomap_\d{3}\.png$", file)
            return json(Dict("ok"=>false, "error"=>"invalid filename"))
        end

        path = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond, "figures", file)
        if !isfile(path)
            return json(Dict("ok"=>false, "error"=>"not found"))
        end

        b64 = Base64.base64encode(read(path))
        json(Dict("ok"=>true, "src"=>"data:image/png;base64,$(b64)"))
    end

    # ─── API: Features de clasificación ICA ──────────────────
    route("/api/ica_features") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        _serve_csv(joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond,
                            "ica_component_features.csv"))
    end

    # ─── API: Fase 6 — Segmentación ───────────────────────────
    route("/api/phase6_segmentation") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        summ_path = joinpath(res_base, "segmentation_summary.json")
        segs_path = joinpath(res_base, "segments_table.csv")
        cov_path  = joinpath(res_base, "channel_coverage.csv")
        log_path  = joinpath(res_base, "pipeline_log.txt")

        # Valores por defecto desde config
        seg_cfg = cfg.segmentation
        ar_cfg  = cfg.artifact_rejection
        epoch_s  = Float64(get(seg_cfg, "epoch_length_s", 2.0))
        overlap_f= Float64(get(seg_cfg, "epoch_overlap",  0.0))
        overlap_s= round(epoch_s * overlap_f, digits=3)
        step_s   = round(epoch_s * (1.0 - overlap_f), digits=3)
        fs_cfg   = Float64(get(cfg.recording, "fs", 500.0))
        samp_ep  = round(Int, epoch_s * fs_cfg)

        default_summary = Dict(
            "n_total"              => 0,
            "n_valid"              => 0,
            "n_rejected"           => 0,
            "retention_pct"        => 0.0,
            "epoch_length_s"       => epoch_s,
            "overlap_s"            => overlap_s,
            "overlap_pct"          => round(overlap_f * 100, digits=1),
            "step_s"               => step_s,
            "n_channels"           => 0,
            "fs"                   => fs_cfg,
            "samples_per_epoch"    => samp_ep,
            "signal_duration_s"    => 0.0,
            "signal_input"         => "filtrada",
            "amp_threshold_uv"     => Float64(get(ar_cfg, "amplitude_threshold_uv", 100.0)),
            "grad_threshold_uv"    => Float64(get(ar_cfg, "gradient_threshold_uv",  50.0)),
            "baseline_method"      => "mean",
            "n_rejected_amplitude" => 0,
            "n_rejected_gradient"  => 0,
            "quality_mean"         => 0.0,
            "quality_median"       => 0.0,
            "quality_threshold"    => 0.5,
            "quality_histogram"    => [],
            "coverage_mean_pct"    => 0.0,
            "coverage_min_pct"     => 0.0,
            "coverage_max_pct"     => 0.0,
            "timestamp"            => "",
            "duration_s"           => 0.0,
        )

        if !isfile(summ_path)
            return json(Dict(
                "ok"          => true,
                "seg_run"     => false,
                "summary"     => default_summary,
                "segments"    => Dict{String,Any}[],
                "channel_coverage" => Dict{String,Any}[],
                "phase_timing"=> Dict("start"=>"","end"=>"","duration"=>""),
            ))
        end

        # ── Parsear segmentation_summary.json ──────────────────
        summary = Dict{String,Any}(default_summary)
        try
            txt = read(summ_path, String)
            for (key, T) in [
                ("n_total",Int), ("n_valid",Int), ("n_rejected",Int),
                ("n_channels",Int), ("samples_per_epoch",Int),
                ("n_rejected_amplitude",Int), ("n_rejected_gradient",Int),
            ]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                m !== nothing && (summary[key] = parse(Int, m.captures[1]))
            end
            for (key, T) in [
                ("retention_pct",Float64), ("epoch_length_s",Float64),
                ("overlap_s",Float64), ("overlap_pct",Float64), ("step_s",Float64),
                ("fs",Float64), ("signal_duration_s",Float64),
                ("amp_threshold_uv",Float64), ("grad_threshold_uv",Float64),
                ("quality_mean",Float64), ("quality_median",Float64),
                ("quality_threshold",Float64), ("coverage_mean_pct",Float64),
                ("coverage_min_pct",Float64), ("coverage_max_pct",Float64),
                ("duration_s",Float64),
            ]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
            end
            for key in ["signal_input","baseline_method","timestamp"]
                m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                m !== nothing && (summary[key] = String(m.captures[1]))
            end
            # histograma: array de [bin, count]
            mh = match(r"\"quality_histogram\"\s*:\s*(\[[^\]]*\])", txt)
            if mh !== nothing
                raw_hist = mh.captures[1]
                pairs = collect(eachmatch(r"\[([0-9.]+),([0-9]+)\]", raw_hist))
                summary["quality_histogram"] = [
                    Dict("bin"=>parse(Float64,p.captures[1]),
                         "count"=>parse(Int,p.captures[2])) for p in pairs
                ]
            end
        catch e
            @warn "segmentation_summary.json parse error: $e"
        end

        # ── Parsear segments_table.csv ──────────────────────────
        segments = Dict{String,Any}[]
        if isfile(segs_path)
            try
                df = CSV.read(segs_path, DataFrame)
                for row in eachrow(df)
                    push!(segments, Dict{String,Any}(
                        "epoch"            => Int(row.epoch),
                        "start_s"          => Float64(row.start_s),
                        "end_s"            => Float64(row.end_s),
                        "duration_s"       => Float64(row.duration_s),
                        "quality"          => round(Float64(row.quality), digits=3),
                        "status"           => string(row.status),
                        "rejection_reason" => string(get(row, :rejection_reason, "")),
                        "max_amp_uv"       => round(Float64(get(row, :max_amp_uv, 0.0)), digits=1),
                        "max_grad_uv"      => round(Float64(get(row, :max_grad_uv, 0.0)), digits=1),
                    ))
                end
            catch e
                @warn "segments_table.csv parse error: $e"
            end
        end

        # ── Parsear channel_coverage.csv ────────────────────────
        ch_coverage = Dict{String,Any}[]
        if isfile(cov_path)
            try
                df = CSV.read(cov_path, DataFrame)
                for row in eachrow(df)
                    push!(ch_coverage, Dict{String,Any}(
                        "channel"      => string(row.channel),
                        "coverage_pct" => Float64(row.coverage_pct),
                    ))
                end
            catch e
                @warn "channel_coverage.csv parse error: $e"
            end
        end

        # ── Timing desde log ────────────────────────────────────
        t_start_str = ""; t_end_str = ""; t_dur = ""
        if isfile(log_path)
            try
                txt = read(log_path, String)
                lines = split(txt, '\n')
                for line in lines
                    if occursin("[5/8] Segmentación", line) || occursin("[5/8] Segment", line)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_start_str = String(m.captures[1]))
                    end
                    if occursin("Duración segmentación", line)
                        m  = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m2 = match(r"(\d+\.?\d*)\s*s", line)
                        m !== nothing && (t_end_str = String(m.captures[1]))
                        m2 !== nothing && (t_dur = String(m2.captures[1]) * " s")
                    end
                end
            catch; end
        end

        json(Dict(
            "ok"             => true,
            "seg_run"        => true,
            "summary"        => summary,
            "segments"       => segments,
            "channel_coverage" => ch_coverage,
            "phase_timing"   => Dict(
                "start"    => t_start_str,
                "end"      => t_end_str,
                "duration" => t_dur,
            ),
        ))
    end

    # ─── API: Fase 7 — Rechazo de Artefactos ──────────────────
    route("/api/phase7_ar") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        summ_path = joinpath(res_base, "artifact_rejection_summary.json")
        rej_path  = joinpath(res_base, "rejected_segments.csv")
        ch_path   = joinpath(res_base, "channel_artifact_summary.csv")
        log_path  = joinpath(res_base, "pipeline_log.txt")

        # Valores por defecto desde config
        ar_cfg  = cfg.artifact_rejection
        seg_cfg = cfg.segmentation
        ar_prof_def = String(get(ar_cfg, "profile", "default"))
        is_ej_def   = ar_prof_def == "eeg_julia"

        default_summary = Dict(
            "profile"                => ar_prof_def,
            "n_total"                => 0,
            "n_valid"                => 0,
            "n_rejected"             => 0,
            "retention_pct"          => 0.0,
            "n_rejected_amplitude"   => 0,
            "n_rejected_gradient"    => 0,
            "min_amplitude_uv"       => Float64(get(ar_cfg, "min_amplitude_uv",        is_ej_def ? -70.0 : -100.0)),
            "max_amplitude_uv"       => Float64(get(ar_cfg, "max_amplitude_uv",        is_ej_def ?  70.0 :  100.0)),
            "amp_threshold_uv"       => Float64(get(ar_cfg, "amplitude_threshold_uv",  100.0)),
            "grad_threshold_uv"      => Float64(get(ar_cfg, "gradient_threshold_uv",    50.0)),
            "use_gradient"           => Bool(get(ar_cfg,   "use_gradient",             !is_ej_def)),
            "n_channels_used"        => Int(get(ar_cfg,    "n_channels_used",           30)),
            "n_channels_total"       => 0,
            "before_event_ms"        => Int(get(ar_cfg,    "before_event_ms",           200)),
            "after_event_ms"         => Int(get(ar_cfg,    "after_event_ms",            300)),
            "before_after_applied"   => false,
            "p2p_mean_uv"            => 0.0,
            "p2p_std_uv"             => 0.0,
            "p2p_max_uv"             => 0.0,
            "p2p_thresh_2sd"         => 0.0,
            "p2p_histogram"          => [],
            "timestamp"              => "",
            "duration_s"             => 0.0,
        )

        if !isfile(summ_path)
            return json(Dict(
                "ok"               => true,
                "ar_run"           => false,
                "summary"          => default_summary,
                "rejected_segments"=> Dict{String,Any}[],
                "channel_summary"  => Dict{String,Any}[],
                "phase_timing"     => Dict("start"=>"","end"=>"","duration"=>""),
            ))
        end

        # ── Parsear artifact_rejection_summary.json ──────────────────────────
        summary = Dict{String,Any}(default_summary)
        try
            txt = read(summ_path, String)
            for key in ["n_total","n_valid","n_rejected",
                        "n_rejected_amplitude","n_rejected_gradient",
                        "n_channels_used","n_channels_total",
                        "before_event_ms","after_event_ms"]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                m !== nothing && (summary[key] = parse(Int, m.captures[1]))
            end
            for key in ["retention_pct","amp_threshold_uv","grad_threshold_uv",
                        "min_amplitude_uv","max_amplitude_uv",
                        "p2p_mean_uv","p2p_std_uv","p2p_max_uv","p2p_thresh_2sd",
                        "duration_s"]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
            end
            for key in ["timestamp","profile"]
                m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                m !== nothing && (summary[key] = String(m.captures[1]))
            end
            # use_gradient / before_after_applied (boolean)
            for key in ["use_gradient","before_after_applied"]
                m = match(Regex("\"$(key)\"\\s*:\\s*(true|false)"), txt)
                m !== nothing && (summary[key] = m.captures[1] == "true")
            end
            # histograma P2P: array de [bin, count]
            mh = match(r"\"p2p_histogram\"\s*:\s*(\[[^\]]*\])", txt)
            if mh !== nothing
                pairs = collect(eachmatch(r"\[([0-9.]+),([0-9]+)\]", mh.captures[1]))
                summary["p2p_histogram"] = [
                    Dict("bin"=>parse(Float64,p.captures[1]),
                         "count"=>parse(Int,p.captures[2])) for p in pairs
                ]
            end
        catch e
            @warn "artifact_rejection_summary.json parse error: $e"
        end

        # ── Labels derivados para el frontend ─────────────────────────────────
        is_ej    = string(get(summary, "profile", "default")) == "eeg_julia"
        use_grad = Bool(get(summary, "use_gradient", !is_ej))
        min_uv   = Float64(get(summary, "min_amplitude_uv", is_ej ? -70.0 : -100.0))
        max_uv   = Float64(get(summary, "max_amplitude_uv", is_ej ?  70.0 :  100.0))
        n_ch_u   = Int(get(summary, "n_channels_used", 30))
        n_ch_t   = Int(get(summary, "n_channels_total", 0))
        ch_desc  = is_ej ? "primeros $(n_ch_u) canales" :
                           (n_ch_t > 0 ? "todos ($(n_ch_t))" : "todos los canales")
        summary["detector_label"] = is_ej ?
            "Amplitud ±$(Int(round(max_uv))) µV" :
            "Amplitud $(Int(round(max_uv))) µV" * (use_grad ? " + Gradiente" : "")
        summary["gradient_label"] = use_grad ?
            "$(Int(round(Float64(get(summary, "grad_threshold_uv", 50.0))))) µV/muestra" :
            "No aplicado"
        summary["channels_desc"] = ch_desc

        # ── Parsear rejected_segments.csv ────────────────────────────────────
        rejected_segs = Dict{String,Any}[]
        if isfile(rej_path)
            try
                df = CSV.read(rej_path, DataFrame)
                for row in eachrow(df)
                    d = Dict{String,Any}(
                        "epoch"               => Int(get(row, :epoch,    0)),
                        "start_s"             => Float64(get(row, :start_s, 0.0)),
                        "end_s"               => Float64(get(row, :end_s,   0.0)),
                        "status"              => string(get(row, :status,   "rejected")),
                        "rejection_reason"    => string(get(row, :rejection_reason, "")),
                        "max_amp_uv"          => round(Float64(get(row, :max_amp_uv,  0.0)), digits=1),
                        "min_amp_uv"          => round(Float64(get(row, :min_amp_uv,  0.0)), digits=1),
                        "max_grad_uv"         => round(Float64(get(row, :max_grad_uv, 0.0)), digits=1),
                        "p2p_uv"              => round(Float64(get(row, :p2p_uv,      0.0)), digits=1),
                        "worst_channel"       => string(get(row, :worst_channel, "")),
                        "channels_violating"  => string(get(row, :channels_violating, "")),
                        "quality"             => round(Float64(get(row, :quality, 0.0)), digits=3),
                    )
                    push!(rejected_segs, d)
                end
            catch e
                @warn "rejected_segments.csv parse error: $e"
            end
        end

        # ── Parsear channel_artifact_summary.csv ─────────────────────────────
        ch_summary = Dict{String,Any}[]
        if isfile(ch_path)
            try
                df = CSV.read(ch_path, DataFrame)
                for row in eachrow(df)
                    push!(ch_summary, Dict{String,Any}(
                        "channel"     => string(get(row, :channel,     "")),
                        "n_bad"       => Int(get(row, :n_bad,      0)),
                        "pct_bad"     => Float64(get(row, :pct_bad,    0.0)),
                        "main_reason" => string(get(row, :main_reason, "")),
                    ))
                end
            catch e
                @warn "channel_artifact_summary.csv parse error: $e"
            end
        end

        # ── Timing desde log ──────────────────────────────────────────────────
        t_start_str = ""; t_end_str = ""; t_dur = ""
        if isfile(log_path)
            try
                txt = read(log_path, String)
                for line in split(txt, '\n')
                    if occursin("[5/8] Segmentación", line) || occursin("[5/8] Segment", line)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_start_str = String(m.captures[1]))
                    end
                    if occursin("Duración segmentación", line)
                        m  = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m2 = match(r"(\d+\.?\d*)\s*s", line)
                        m !== nothing && (t_end_str = String(m.captures[1]))
                        m2 !== nothing && (t_dur = String(m2.captures[1]) * " s")
                    end
                end
            catch; end
        end

        json(Dict(
            "ok"                => true,
            "ar_run"            => true,
            "summary"           => summary,
            "rejected_segments" => rejected_segs,
            "channel_summary"   => ch_summary,
            "phase_timing"      => Dict(
                "start"    => t_start_str,
                "end"      => t_end_str,
                "duration" => t_dur,
            ),
        ))
    end

    # ─── API: Fase 7 — Señal de un epoch concreto ─────────────
    route("/api/phase7_epoch_signal") do
        subj   = string(get(getpayload(), :subj,    "M05"))
        sess   = string(get(getpayload(), :sess,    "T2"))
        cond   = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        ep_str = string(get(getpayload(), :epoch,   "1"))
        ch_req = string(get(getpayload(), :channel, ""))

        epoch_idx = tryparse(Int, ep_str)
        epoch_idx === nothing && return json(Dict("ok"=>false,"error"=>"invalid epoch"))

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        sig_path  = joinpath(res_base, "ica_signal_after.csv")
        summ_path = joinpath(res_base, "segmentation_summary.json")

        isfile(sig_path) || return json(Dict("ok"=>false,"error"=>"no signal file"))

        try
            df = CSV.read(sig_path, DataFrame)
            # Epoch length from segmentation summary (default 1.0 s)
            epoch_s = 1.0
            if isfile(summ_path)
                txt2 = read(summ_path, String)
                m2 = match(r"\"epoch_length_s\"\s*:\s*([0-9.]+)", txt2)
                m2 !== nothing && (epoch_s = parse(Float64, m2.captures[1]))
            end
            t_start = (epoch_idx - 1) * epoch_s
            t_end   = epoch_idx       * epoch_s

            t_col   = Float64.(df.t_s)
            mask    = (t_col .>= t_start .- 1e-4) .& (t_col .< t_end .+ 1e-4)
            sub_df  = df[mask, :]
            isempty(sub_df) && return json(Dict("ok"=>false,"error"=>"epoch out of range"))

            ch_names_all = [string(c) for c in names(df) if string(c) != "t_s"]
            isempty(ch_names_all) && return json(Dict("ok"=>false,"error"=>"no channels"))

            # Resolve requested channel
            ch_idx = 1
            if !isempty(ch_req)
                idx2 = findfirst(==(ch_req), ch_names_all)
                idx2 !== nothing && (ch_idx = idx2)
            end
            ch_name = ch_names_all[ch_idx]
            ch_sym  = Symbol(ch_name)
            values  = hasproperty(sub_df, ch_sym) ?
                      round.(Float64.(getproperty(sub_df, ch_sym)), digits=3) : Float64[]

            return json(Dict(
                "ok"       => true,
                "epoch"    => epoch_idx,
                "channel"  => ch_name,
                "t_start"  => t_start,
                "t_end"    => t_end,
                "epoch_s"  => epoch_s,
                "times"    => round.(Float64.(sub_df.t_s), digits=4),
                "values"   => values,
                "ch_names" => ch_names_all,
            ))
        catch e
            @warn "phase7_epoch_signal error: $e"
            return json(Dict("ok"=>false,"error"=>string(e)))
        end
    end

    # ─── API: Fase 8 — Análisis Espectral ────────────────────
    route("/api/phase8_spectral") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        res_base  = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        summ_path = joinpath(res_base, "spectral_summary.json")
        bp_path   = joinpath(res_base, "band_power_summary.csv")
        psd_path  = joinpath(res_base, "psd_by_channel.csv")
        reg_path  = joinpath(res_base, "regional_psd.csv")
        idx_path  = joinpath(res_base, "spectral_indices.csv")
        log_path  = joinpath(res_base, "pipeline_log.txt")

        spectral_run = isfile(bp_path)
        if !spectral_run
            return json(Dict(
                "ok"               => true,
                "spectral_run"     => false,
                "summary"          => Dict{String,Any}(),
                "band_power"       => Dict{String,Any}[],
                "regional_psd"     => Dict{String,Any}[],
                "spectral_indices" => Dict{String,Any}[],
                "psd_global"       => Dict("freqs"=>Float64[],"power"=>Float64[]),
                "phase_timing"     => Dict("start"=>"","end"=>"","duration"=>""),
            ))
        end

        # ── Defaults from config ──────────────────────────────────────────────
        summary = Dict{String,Any}(
            "method"         => "fft_hamming_taper",
            "window"         => "hamming_taper",
            "window_pct"     => Float64(get(cfg.spectral, "window_pct", 10.0)),
            "nfft"           => Int(get(cfg.spectral, "nfft", 512)),
            "epoch_length_s" => Float64(get(cfg.segmentation, "epoch_length_s",
                                    get(cfg.segmentation, "segment_length_seconds", 1.0))),
            "n_epochs"       => 0,
            "fs"             => Float64(get(cfg.recording, "fs",
                                    get(cfg.recording, "sampling_rate", 500.0))),
            "delta_f"        => 0.0,
            "freq_range_lo"  => 0.0,
            "freq_range_hi"  => 250.0,
            "n_freq_bins"    => 0,
            "n_channels"     => 0,
            "reference"      => String(get(cfg.recording, "reference", "average")),
            "profile"        => String(get(cfg.filtering, "profile", "default")),
            "total_power_uv2"=> 0.0,
            "timestamp"      => "",
        )

        # ── Parse spectral_summary.json ───────────────────────────────────────
        if isfile(summ_path)
            try
                txt = read(summ_path, String)
                for key in ["nfft","n_epochs","n_freq_bins","n_channels"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                    m !== nothing && (summary[key] = parse(Int, m.captures[1]))
                end
                for key in ["window_pct","epoch_length_s","fs","delta_f",
                            "freq_range_lo","freq_range_hi","total_power_uv2"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                    m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                end
                for key in ["method","window","reference","profile","timestamp"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                    m !== nothing && (summary[key] = String(m.captures[1]))
                end
                # Flat band stats: {BAND}_{stat}
                for bname in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    for stat in ["mean","std","median","min","max","pct"]
                        key = "$(bname)_$(stat)"
                        m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                        m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                    end
                end
            catch e
                @warn "spectral_summary.json parse error: $e"
            end
        end

        # ── Parse band_power_summary.csv ──────────────────────────────────────
        band_power = Dict{String,Any}[]
        try
            df = CSV.read(bp_path, DataFrame)
            for row in eachrow(df)
                d = Dict{String,Any}("channel" => string(row.channel))
                for col in names(df)
                    col == "channel" && continue
                    v = getproperty(row, Symbol(col))
                    d[col] = v isa AbstractFloat ? round(Float64(v), digits=6) : v
                end
                push!(band_power, d)
            end
        catch e
            @warn "band_power_summary.csv parse error: $e"
        end

        # ── Compute global mean PSD (0–80 Hz) from psd_by_channel.csv ────────
        psd_global = Dict("freqs"=>Float64[], "power"=>Float64[])
        if isfile(psd_path)
            try
                df = CSV.read(psd_path, DataFrame)
                freq_groups = Dict{Float64, Vector{Float64}}()
                for row in eachrow(df)
                    f = Float64(row.freq_hz)
                    f > 80.0 && continue     # Limit display range
                    p = Float64(row.power_uv2)
                    v = get!(freq_groups, f, Float64[])
                    push!(v, p)
                end
                freq_sorted = sort(collect(keys(freq_groups)))
                psd_global = Dict(
                    "freqs" => [round(f, digits=3) for f in freq_sorted],
                    "power" => [round(mean(freq_groups[f]), digits=9) for f in freq_sorted],
                )
            catch e
                @warn "psd_by_channel.csv global parse error: $e"
            end
        end

        # ── Parse regional_psd.csv ────────────────────────────────────────────
        regional_psd = Dict{String,Any}[]
        if isfile(reg_path)
            try
                df = CSV.read(reg_path, DataFrame)
                for row in eachrow(df)
                    push!(regional_psd, Dict{String,Any}(
                        "region"     => string(row.region),
                        "band"       => string(row.band),
                        "mean_power" => round(Float64(row.mean_power), digits=6),
                        "std_power"  => round(Float64(row.std_power),  digits=6),
                        "n_channels" => Int(row.n_channels),
                    ))
                end
            catch e
                @warn "regional_psd.csv parse error: $e"
            end
        end

        # ── Parse spectral_indices.csv ────────────────────────────────────────
        spectral_indices = Dict{String,Any}[]
        if isfile(idx_path)
            try
                df = CSV.read(idx_path, DataFrame)
                for row in eachrow(df)
                    push!(spectral_indices, Dict{String,Any}(
                        "channel"        => string(row.channel),
                        "alpha_theta"    => round(Float64(get(row, :alpha_theta,    0.0)), digits=3),
                        "beta_alpha"     => round(Float64(get(row, :beta_alpha,     0.0)), digits=3),
                        "theta_beta"     => round(Float64(get(row, :theta_beta,     0.0)), digits=3),
                        "gamma_alpha"    => round(Float64(get(row, :gamma_alpha,    0.0)), digits=3),
                        "peak_alpha_hz"  => round(Float64(get(row, :peak_alpha_hz,  0.0)), digits=2),
                        "peak_alpha_uv2" => round(Float64(get(row, :peak_alpha_uv2, 0.0)), digits=6),
                    ))
                end
            catch e
                @warn "spectral_indices.csv parse error: $e"
            end
        end

        # ── Timing from log ───────────────────────────────────────────────────
        t_start = ""; t_end = ""; t_dur = ""
        if isfile(log_path)
            try
                txt = read(log_path, String)
                for line in split(txt, '\n')
                    if (occursin("[6/8]", line) || occursin("Espectral", line)) && isempty(t_start)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_start = String(m.captures[1]))
                    end
                    if occursin("[7/8]", line) && isempty(t_end)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_end = String(m.captures[1]))
                    end
                end
            catch; end
        end

        json(Dict(
            "ok"               => true,
            "spectral_run"     => true,
            "summary"          => summary,
            "band_power"       => band_power,
            "regional_psd"     => regional_psd,
            "spectral_indices" => spectral_indices,
            "psd_global"       => psd_global,
            "phase_timing"     => Dict("start"=>t_start,"end"=>t_end,"duration"=>t_dur),
        ))
    end

    # ─── API: Fase 9 — Conectividad wPLI ─────────────────────
    route("/api/phase9_wpli") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        band = string(get(getpayload(), :band, "ALPHA"))

        res_base   = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        edges_path = joinpath(res_base, "connectivity_edges.csv")
        summ_path  = joinpath(res_base, "connectivity_summary.json")
        nm_path    = joinpath(res_base, "network_metrics.csv")
        log_path   = joinpath(res_base, "pipeline_log.txt")
        mat_path   = joinpath(res_base, "wpli_$(band).csv")

        conn_run = isfile(edges_path) || isfile(mat_path)
        if !conn_run
            return json(Dict(
                "ok"              => true,
                "conn_run"        => false,
                "summary"         => Dict{String,Any}(),
                "matrix"          => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"band"=>band),
                "edges"           => Dict{String,Any}[],
                "network_metrics" => Dict{String,Any}[],
                "phase_timing"    => Dict("start"=>"","end"=>"","duration"=>""),
            ))
        end

        # ── Parse connectivity_summary.json ────────────────────────────────────
        summary = Dict{String,Any}(
            "method"        => "wpli",
            "space"         => "sensor",
            "n_channels"    => 0,
            "n_epochs_used" => 0,
            "threshold"     => 0.1,
            "timestamp"     => "",
        )
        if isfile(summ_path)
            try
                txt = read(summ_path, String)
                for key in ["n_channels","n_epochs_used"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                    m !== nothing && (summary[key] = parse(Int, m.captures[1]))
                end
                for key in ["threshold"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                    m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                end
                for key in ["method","space","timestamp"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                    m !== nothing && (summary[key] = String(m.captures[1]))
                end
                # Per-band stats: BAND_stat
                for bname in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    for stat in ["mean","std","median","max","density"]
                        key = "$(bname)_$(stat)"
                        m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                        m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                    end
                    for stat in ["n_edges","n_above"]
                        key = "$(bname)_$(stat)"
                        m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                        m !== nothing && (summary[key] = parse(Int, m.captures[1]))
                    end
                end
            catch e
                @warn "connectivity_summary.json parse error: $e"
            end
        end

        # ── Parse wPLI matrix for selected band ────────────────────────────────
        matrix = Dict{String,Any}("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"band"=>band)
        if isfile(mat_path)
            try
                df  = CSV.read(mat_path, DataFrame)
                chs = [string(row.channel) for row in eachrow(df)]
                val_cols = [c for c in names(df) if c != "channel"]
                vals = Vector{Vector{Float64}}()
                for row in eachrow(df)
                    push!(vals, [round(Float64(getproperty(row, Symbol(c))), digits=4) for c in val_cols])
                end
                matrix = Dict{String,Any}("channels"=>chs,"values"=>vals,"band"=>band)
            catch e
                @warn "wpli matrix parse error band=$(band): $e"
            end
        end

        # ── Parse connectivity_edges.csv — top 30 for selected band ───────────
        edges = Dict{String,Any}[]
        if isfile(edges_path)
            try
                df      = CSV.read(edges_path, DataFrame)
                df_band = filter(row -> string(row.band) == band, df)
                sort!(df_band, :wpli; rev=true)
                n_top   = min(30, nrow(df_band))
                for row in eachrow(df_band[1:n_top, :])
                    wv = Float64(row.wpli)
                    push!(edges, Dict{String,Any}(
                        "ch_a"     => string(row.ch_a),
                        "ch_b"     => string(row.ch_b),
                        "band"     => string(row.band),
                        "wpli"     => round(wv, digits=4),
                        "rank"     => Int(row.rank),
                        "fisher_z" => round(atanh(min(wv, 0.9999)), digits=4),
                    ))
                end
            catch e
                @warn "connectivity_edges.csv parse error: $e"
            end
        end

        # ── Parse network_metrics.csv ─────────────────────────────────────────
        network_metrics = Dict{String,Any}[]
        if isfile(nm_path)
            try
                df = CSV.read(nm_path, DataFrame)
                for row in eachrow(df)
                    push!(network_metrics, Dict{String,Any}(
                        "channel"       => string(row.channel),
                        "strength"      => round(Float64(row.strength),      digits=4),
                        "degree"        => Int(row.degree),
                        "norm_strength" => round(Float64(row.norm_strength), digits=4),
                    ))
                end
            catch e
                @warn "network_metrics.csv parse error: $e"
            end
        end

        # ── Timing from log ───────────────────────────────────────────────────
        t_start = ""; t_end = ""; t_dur = ""
        if isfile(log_path)
            try
                txt = read(log_path, String)
                for line in split(txt, '\n')
                    if (occursin("[7/8]", line) || occursin("wPLI", line) ||
                        occursin("Conectividad", line)) && isempty(t_start)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_start = String(m.captures[1]))
                    end
                    if occursin("[8/8]", line) && isempty(t_end)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_end = String(m.captures[1]))
                    end
                end
            catch; end
        end

        json(Dict(
            "ok"              => true,
            "conn_run"        => true,
            "summary"         => summary,
            "matrix"          => matrix,
            "edges"           => edges,
            "network_metrics" => network_metrics,
            "phase_timing"    => Dict("start"=>t_start,"end"=>t_end,"duration"=>t_dur),
        ))
    end

    # ─── API: Fase 10 — Surrogates / Inferencia ──────────────
    route("/api/phase10_surrogates") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))
        band = string(get(getpayload(), :band, "ALPHA"))

        res_base   = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        summ_path  = joinpath(res_base, "surrogate_summary.json")
        sig_path   = joinpath(res_base, "significant_connections.csv")
        qc_path    = joinpath(res_base, "surrogate_quality.csv")
        null_path  = joinpath(res_base, "surrogate_null_stats_$(band).csv")
        obs_path   = joinpath(res_base, "wpli_observed_$(band).csv")
        pval_path  = joinpath(res_base, "wpli_pvalues_$(band).csv")
        log_path   = joinpath(res_base, "pipeline_log.txt")

        surr_run = isfile(summ_path) || isfile(sig_path)
        if !surr_run
            return json(Dict(
                "ok"               => true,
                "surr_run"         => false,
                "summary"          => Dict{String,Any}(),
                "matrix_obs"       => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"band"=>band),
                "sig_connections"  => Dict{String,Any}[],
                "null_stats"       => Dict{String,Any}[],
                "surrogate_quality"=> Dict{String,Any}[],
                "pvalue_hist"      => Dict{String,Any}[],
                "phase_timing"     => Dict("start"=>"","end"=>"","duration"=>""),
            ))
        end

        # ── Parse surrogate_summary.json ──────────────────────────────────────
        summary = Dict{String,Any}(
            "method"        => "phase_shuffle",
            "n_surrogates"  => 200,
            "alpha"         => 0.05,
            "fdr_method"    => "bh",
            "seed"          => 42,
            "n_channels"    => 0,
            "n_total_pairs" => 0,
            "n_sig_total"   => 0,
            "timestamp"     => "",
        )
        if isfile(summ_path)
            try
                txt = read(summ_path, String)
                for key in ["n_surrogates","seed","n_channels","n_total_pairs","n_sig_total"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                    m !== nothing && (summary[key] = parse(Int, m.captures[1]))
                end
                for key in ["alpha"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                    m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                end
                for key in ["method","fdr_method","timestamp"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                    m !== nothing && (summary[key] = String(m.captures[1]))
                end
                # Per-band: BAND_n_sig, BAND_pct_sig, BAND_mean_p, BAND_fdr_thr
                for bname in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    for stat in ["pct_sig","mean_p","fdr_thr"]
                        key = "$(bname)_$(stat)"
                        m = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                        m !== nothing && (summary[key] = parse(Float64, m.captures[1]))
                    end
                    for stat in ["n_sig"]
                        key = "$(bname)_$(stat)"
                        m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                        m !== nothing && (summary[key] = parse(Int, m.captures[1]))
                    end
                end
            catch e
                @warn "surrogate_summary.json parse error: $e"
            end
        end

        # ── Parse significant_connections.csv ─────────────────────────────────
        sig_connections = Dict{String,Any}[]
        if isfile(sig_path)
            try
                df      = CSV.read(sig_path, DataFrame)
                df_band = filter(r -> string(r.band) == band, df)
                sort!(df_band, :q_value)
                for row in eachrow(df_band[1:min(50,nrow(df_band)), :])
                    push!(sig_connections, Dict{String,Any}(
                        "ch_a"     => string(row.ch_a),
                        "ch_b"     => string(row.ch_b),
                        "band"     => string(row.band),
                        "wpli_obs" => round(Float64(row.wpli_obs), digits=4),
                        "p_value"  => round(Float64(row.p_value),  digits=4),
                        "q_value"  => round(Float64(row.q_value),  digits=4),
                        "z_score"  => round(Float64(row.z_score),  digits=3),
                    ))
                end
            catch e
                @warn "significant_connections.csv parse error: $e"
            end
        end

        # ── Parse wpli_observed_{band}.csv ────────────────────────────────────
        matrix_obs = Dict{String,Any}("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"band"=>band)
        if isfile(obs_path)
            try
                df  = CSV.read(obs_path, DataFrame)
                chs = [string(row.channel) for row in eachrow(df)]
                val_cols = [c for c in names(df) if c != "channel"]
                vals = Vector{Vector{Float64}}()
                for row in eachrow(df)
                    push!(vals, [round(Float64(getproperty(row, Symbol(c))), digits=4) for c in val_cols])
                end
                matrix_obs = Dict{String,Any}("channels"=>chs,"values"=>vals,"band"=>band)
            catch e
                @warn "wpli_observed parse error for band=$(band): $e"
            end
        end

        # ── Parse surrogate_null_stats_{band}.csv ─────────────────────────────
        null_stats = Dict{String,Any}[]
        if isfile(null_path)
            try
                df = CSV.read(null_path, DataFrame)
                for row in eachrow(df)
                    push!(null_stats, Dict{String,Any}(
                        "ch_a"      => string(row.ch_a),
                        "ch_b"      => string(row.ch_b),
                        "wpli_obs"  => round(Float64(row.wpli_obs),   digits=4),
                        "null_mean" => round(Float64(row.null_mean), digits=4),
                        "null_std"  => round(Float64(row.null_std),  digits=4),
                        "p_value"   => round(Float64(row.p_value),   digits=4),
                        "q_value"   => round(Float64(row.q_value),   digits=4),
                        "z_score"   => round(Float64(row.z_score),   digits=3),
                    ))
                end
            catch e
                @warn "surrogate_null_stats parse error: $e"
            end
        end

        # ── Compute p-value histogram from wpli_pvalues_{band}.csv ───────────
        pvalue_hist = Dict{String,Any}[]
        if isfile(pval_path)
            try
                df       = CSV.read(pval_path, DataFrame)
                val_cols = [c for c in names(df) if c != "channel"]
                all_p    = Float64[]
                for (i, row) in enumerate(eachrow(df))
                    for (j, col) in enumerate(val_cols)
                        j <= i && continue   # upper triangle only
                        v = getproperty(row, Symbol(col))
                        push!(all_p, Float64(v))
                    end
                end
                n_bins = 20
                bins   = zeros(Int, n_bins)
                for p in all_p
                    b = min(n_bins, max(1, ceil(Int, p * n_bins)))
                    bins[b] += 1
                end
                for i in 1:n_bins
                    push!(pvalue_hist, Dict{String,Any}(
                        "lo" => round((i-1)/n_bins, digits=3),
                        "hi" => round(i/n_bins, digits=3),
                        "count" => bins[i],
                    ))
                end
            catch e
                @warn "p-value histogram error: $e"
            end
        end

        # ── Parse surrogate_quality.csv ───────────────────────────────────────
        surrogate_quality = Dict{String,Any}[]
        if isfile(qc_path)
            try
                df = CSV.read(qc_path, DataFrame)
                for row in eachrow(df)
                    push!(surrogate_quality, Dict{String,Any}(
                        "band"           => string(row.band),
                        "n_surrogates"   => Int(row.n_surrogates),
                        "n_sig"          => Int(row.n_sig),
                        "n_total"        => Int(row.n_total),
                        "pct_sig"        => round(Float64(row.pct_sig),        digits=2),
                        "mean_p"         => round(Float64(row.mean_p),         digits=4),
                        "fdr_threshold"  => round(Float64(row.fdr_threshold),  digits=4),
                        "mean_null_mean" => round(Float64(row.mean_null_mean), digits=4),
                        "mean_null_std"  => round(Float64(row.mean_null_std),  digits=4),
                        "obs_mean"       => round(Float64(row.obs_mean),       digits=4),
                        "obs_max"        => round(Float64(row.obs_max),        digits=4),
                    ))
                end
            catch e
                @warn "surrogate_quality.csv parse error: $e"
            end
        end

        # ── Timing from log ───────────────────────────────────────────────────
        t_start = ""; t_end = ""; t_dur = ""
        if isfile(log_path)
            try
                txt = read(log_path, String)
                for line in split(txt, '\n')
                    if occursin("[SUR]", line) && isempty(t_start)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_start = String(m.captures[1]))
                    end
                    if occursin("[8/8]", line) && isempty(t_end)
                        m = match(r"\[(\d{2}:\d{2}:\d{2})\]", line)
                        m !== nothing && (t_end = String(m.captures[1]))
                    end
                end
            catch; end
        end

        json(Dict(
            "ok"               => true,
            "surr_run"         => true,
            "summary"          => summary,
            "matrix_obs"       => matrix_obs,
            "sig_connections"  => sig_connections,
            "null_stats"       => null_stats,
            "surrogate_quality"=> surrogate_quality,
            "pvalue_hist"      => pvalue_hist,
            "phase_timing"     => Dict("start"=>t_start,"end"=>t_end,"duration"=>t_dur),
        ))
    end

    # ─── Legacy API (mantener compatibilidad) ─────────────────

    route("/api/subjects") do
        dirs = isdir(res_root) ?
               filter(d -> isdir(joinpath(res_root, d)) &&
                            d ∉ ["subjects", "transversal", "longitudinal",
                                 "qc", "logs", "deprecated"],
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

    # ─── API: Fase 11 — Resultados Finales ──────────────────────
    route("/api/phase11_summary") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        base = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        fe(f) = isfile(joinpath(base, f))

        run = fe("overview.csv") || fe("pipeline_log.txt")
        if !run
            return json(Dict("ok" => true, "run" => false,
                             "subject" => subj, "session" => sess, "cond" => cond))
        end

        # ── Recording info ──────────────────────────────────────
        recording = Dict{String,Any}(
            "n_channels" => 0, "n_bad_ch" => 0, "bad_channels" => "",
            "duration_s" => 0.0, "fs" => 500.0,
            "n_epochs" => 0, "n_epochs_accepted" => 0, "pct_epochs_retained" => 0.0,
        )
        if fe("overview.csv")
            try
                df = CSV.read(joinpath(base, "overview.csv"), DataFrame)
                r  = df[1, :]
                for (k, sym) in [("n_channels",:n_channels),("n_bad_ch",:n_bad_ch),
                                  ("duration_s",:duration_s),("fs",:fs)]
                    hasproperty(df, sym) || continue
                    v = r[sym]
                    v isa Missing || (recording[k] = Float64(v))
                end
                hasproperty(df, :bad_channels) &&
                    (recording["bad_channels"] = string(r[:bad_channels]))
            catch; end
        end

        # ── Epoch info from log ─────────────────────────────────
        log_txt = fe("pipeline_log.txt") ? read(joinpath(base, "pipeline_log.txt"), String) : ""
        if !isempty(log_txt)
            m = match(r"Epochs aceptados[:\s]+(\d+)/(\d+)", log_txt)
            if m !== nothing
                n_acc = parse(Int, m.captures[1])
                n_tot = parse(Int, m.captures[2])
                recording["n_epochs_accepted"] = n_acc
                recording["n_epochs"]          = n_tot
                recording["pct_epochs_retained"] = n_tot > 0 ?
                    round(100.0 * n_acc / n_tot, digits=1) : 0.0
            end
        end

        # ── ICA info ────────────────────────────────────────────
        ica = Dict{String,Any}("n_comp" => 0, "n_rejected" => 0,
                                "variance_retained" => 1.0, "has_features" => false)
        if fe("ica_summary.json")
            try
                txt = read(joinpath(base, "ica_summary.json"), String)
                for key in ["n_comp","n_rej","n_topomaps"]
                    m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                    m !== nothing &&
                        (ica[key == "n_rej" ? "n_rejected" : key] = parse(Int, m.captures[1]))
                end
                m = match(r"\"variance_retained\"\s*:\s*([0-9.eE+\-]+)", txt)
                m !== nothing && (ica["variance_retained"] = parse(Float64, m.captures[1]))
                m = match(r"\"has_features\"\s*:\s*(true|false)", txt)
                m !== nothing && (ica["has_features"] = m.captures[1] == "true")
            catch; end
        end

        # ── Spectral (band_power_summary.csv) ───────────────────
        band_powers    = Dict{String,Any}[]
        dominant_band  = ""
        dominant_power = -1.0
        if fe("band_power_summary.csv")
            try
                df = CSV.read(joinpath(base, "band_power_summary.csv"), DataFrame)
                for row in eachrow(df)
                    bname = string(row[:band])
                    pval  = hasproperty(df, :mean_power) ? Float64(row[:mean_power]) :
                            hasproperty(df, :power_db)   ? Float64(row[:power_db])   : 0.0
                    rval  = hasproperty(df, :rel_power)  ? Float64(row[:rel_power])  : 0.0
                    push!(band_powers, Dict("band" => bname,
                        "power"     => round(pval, digits=3),
                        "rel_power" => round(rval, digits=4)))
                    if pval > dominant_power
                        dominant_power = pval; dominant_band = bname
                    end
                end
            catch; end
        end

        # ── Connectivity summary ─────────────────────────────────
        connectivity = Dict{String,Any}("best_band" => "", "mean_wpli" => Dict{String,Any}())
        if fe("connectivity_summary.json")
            try
                txt = read(joinpath(base, "connectivity_summary.json"), String)
                m = match(r"\"best_band\"\s*:\s*\"([^\"]+)\"", txt)
                m !== nothing && (connectivity["best_band"] = m.captures[1])
                for bn in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    mk = match(Regex("\"$(bn)_mean_wpli\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                    mk !== nothing &&
                        (connectivity["mean_wpli"][bn] = parse(Float64, mk.captures[1]))
                end
            catch; end
        end

        # ── Surrogate summary ────────────────────────────────────
        surrogates = Dict{String,Any}("n_sig_total" => 0, "best_band" => "",
                                       "has_surrogates" => false)
        if fe("surrogate_summary.json")
            try
                txt = read(joinpath(base, "surrogate_summary.json"), String)
                surrogates["has_surrogates"] = true
                m = match(r"\"n_sig_total\"\s*:\s*([0-9]+)", txt)
                m !== nothing && (surrogates["n_sig_total"] = parse(Int, m.captures[1]))
                m = match(r"\"best_band\"\s*:\s*\"([^\"]+)\"", txt)
                m !== nothing && (surrogates["best_band"] = m.captures[1])
                for bn in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    for (sfx, T) in [("pct_sig",Float64),("n_sig",Int)]
                        key = "$(bn)_$(sfx)"
                        mk = match(Regex("\"$(key)\"\\s*:\\s*([0-9.eE+\\-]+)"), txt)
                        mk !== nothing && (surrogates[key] =
                            T == Int ? parse(Int, mk.captures[1]) :
                                       parse(Float64, mk.captures[1]))
                    end
                end
            catch; end
        end

        # ── Top significant connections ──────────────────────────
        top_connections = Dict{String,Any}[]
        if fe("significant_connections.csv")
            try
                df = CSV.read(joinpath(base, "significant_connections.csv"), DataFrame)
                sort!(df, :q_value)
                for row in eachrow(df[1:min(10, nrow(df)), :])
                    push!(top_connections, Dict{String,Any}(
                        "ch_a"    => string(row.ch_a),
                        "ch_b"    => string(row.ch_b),
                        "band"    => string(row.band),
                        "wpli"    => round(Float64(row.wpli_obs), digits=4),
                        "q_value" => round(Float64(row.q_value),  digits=4),
                    ))
                end
            catch; end
        end

        # ── Pipeline timing ──────────────────────────────────────
        timing = Dict{String,Any}("start" => "", "duration" => "", "phases_done" => 0)
        if !isempty(log_txt)
            m = match(r"NeuroMIND pipeline — (\S+)", log_txt)
            m !== nothing && (timing["start"] = m.captures[1])
            m2 = match(r"completado en ([\d.]+) s", log_txt)
            m2 !== nothing && (timing["duration"] = m2.captures[1] * " s")
        end

        # ── Phase statuses ───────────────────────────────────────
        statuses = Dict{String,String}(
            "0"  => "completed",
            "1"  => fe("overview.csv")           ? "completed" : "pending",
            "2"  => fe("overview.csv")           ? "completed" : "pending",
            "3"  => fe("channel_statistics.csv") ? "completed" : "pending",
            "4"  => !isempty(log_txt)            ? "completed" : "pending",
            "5"  => fe("ica_summary.json")       ? "completed" : "pending",
            "6"  => !isempty(log_txt)            ? "completed" : "pending",
            "7"  => !isempty(log_txt)            ? "completed" : "pending",
            "8"  => fe("psd_by_channel.csv")     ? "completed" : "pending",
            "9"  => fe("wpli_ALPHA.csv") || fe("connectivity_summary.json") ?
                    "completed" : "pending",
            "10" => fe("surrogate_summary.json") ? "completed" : "pending",
        )
        n_done = count(v -> v == "completed", values(statuses))
        timing["phases_done"] = n_done

        # ── Quality score ────────────────────────────────────────
        score_parts = Float64[]
        n_ch  = max(Int(round(get(recording, "n_channels", 31.0))), 1)
        n_bad = Int(round(get(recording, "n_bad_ch", 0.0)))
        push!(score_parts, max(0.0, 1.0 - n_bad / n_ch * 2.0))

        pct_ep = Float64(get(recording, "pct_epochs_retained", 80.0))
        push!(score_parts, clamp(pct_ep / 100.0, 0.0, 1.0))

        var_ret = Float64(get(ica, "variance_retained", 1.0))
        push!(score_parts, clamp(var_ret, 0.0, 1.0))

        push!(score_parts, min(n_done / 10.0, 1.0))

        quality_score = round(mean(score_parts), digits=3)

        json(Dict(
            "ok"              => true,
            "run"             => true,
            "subject"         => subj,
            "session"         => sess,
            "cond"            => cond,
            "recording"       => recording,
            "ica"             => ica,
            "band_powers"     => band_powers,
            "dominant_band"   => dominant_band,
            "connectivity"    => connectivity,
            "surrogates"      => surrogates,
            "top_connections" => top_connections,
            "timing"          => timing,
            "statuses"        => statuses,
            "quality_score"   => quality_score,
        ))
    end

    # ─── API: Fase 12 — Exportación / Reporte ────────────────────
    route("/api/phase12_files") do
        subj = string(get(getpayload(), :subj, "M05"))
        sess = string(get(getpayload(), :sess, "T2"))
        cond = _normalize_cond(string(get(getpayload(), :cond, "EC")))

        base    = joinpath(bids_root, "sub-$(subj)", "ses-$(sess)", cond)
        fig_dir = joinpath(base, "figures")
        run     = isdir(base)

        if !run
            return json(Dict("ok" => true, "run" => false, "categories" => Dict{String,Any}[],
                             "total_files" => 0, "total_size_kb" => 0.0,
                             "subject" => subj, "session" => sess, "cond" => cond))
        end

        BANDS = ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]

        categories_def = [
            ("Registros y QC", [
                "overview.csv","qc_summary.csv","channel_statistics.csv",
                "pipeline_log.txt","config_snapshot.toml"]),
            ("ICA", [
                "ica_summary.json","ica_components.csv","ica_component_features.csv",
                "ica_mixing_matrix.csv","ica_unmixing_matrix.csv",
                "ica_activations.csv","ica_signal_before.csv","ica_signal_after.csv"]),
            ("Espectral", ["psd_by_channel.csv","band_power_summary.csv"]),
            ("Conectividad wPLI", ["connectivity_summary.json","connectivity_edges.csv"]),
            ("Surrogates", [
                "surrogate_summary.json","surrogate_quality.csv",
                "significant_connections.csv"]),
            ("wPLI por banda", vcat([
                ["wpli_$(b).csv","wpli_observed_$(b).csv","wpli_pvalues_$(b).csv",
                 "wpli_qvalues_$(b).csv","wpli_significant_$(b).csv",
                 "surrogate_null_stats_$(b).csv"]
                for b in BANDS]...)),
        ]

        # Figures (from figures/ subdir)
        fig_files = isdir(fig_dir) ? readdir(fig_dir) : String[]
        push!(categories_def, ("Figuras", fig_files))

        result_cats = Dict{String,Any}[]
        total_files = 0
        total_size  = 0.0

        for (cat_name, file_list) in categories_def
            is_fig = cat_name == "Figuras"
            files_info = Dict{String,Any}[]
            for f in file_list
                fpath  = is_fig ? joinpath(fig_dir, f) : joinpath(base, f)
                exists = isfile(fpath)
                sz     = exists ? round(filesize(fpath) / 1024, digits=1) : 0.0
                ext    = endswith(f,".csv") ? "csv" : endswith(f,".json") ? "json" :
                         endswith(f,".png") ? "png" : endswith(f,".txt")  ? "txt"  :
                         endswith(f,".toml") ? "toml" : "other"
                push!(files_info, Dict{String,Any}(
                    "name"    => f,
                    "exists"  => exists,
                    "size_kb" => sz,
                    "type"    => ext,
                    "path"    => is_fig ? "figures/$(f)" : f,
                ))
                if exists; total_files += 1; total_size += sz; end
            end
            n_ok = count(f -> f["exists"], files_info)
            push!(result_cats, Dict{String,Any}(
                "name"      => cat_name,
                "files"     => files_info,
                "n_total"   => length(files_info),
                "n_present" => n_ok,
            ))
        end

        json(Dict(
            "ok"            => true,
            "run"           => true,
            "subject"       => subj,
            "session"       => sess,
            "cond"          => cond,
            "categories"    => result_cats,
            "total_files"   => total_files,
            "total_size_kb" => round(total_size, digits=1),
        ))
    end

    # ─── API: Fase 13 — Evaluación Transversal ──────────────────
    route("/api/phase13_transversal") do
        cond_raw = string(get(getpayload(), :cond, "EC"))
        cond     = _normalize_cond(cond_raw)
        band     = uppercase(string(get(getpayload(), :band, "ALPHA")))

        # El script guarda en .../group/transversal/EC/ o /EO/
        cond_short = cond == "eyesclosed" ? "EC" : (cond == "eyesopen" ? "EO" : uppercase(cond_raw))
        grp_dir    = joinpath(res_root, "transversal", cond_short)

        # Función auxiliar de existencia en grp_dir
        gfe(f) = isfile(joinpath(grp_dir, f))

        empty_resp = Dict{String,Any}(
            "ok"            => true,
            "run"           => false,
            "cond"          => cond_short,
            "band"          => band,
            "summary"       => Dict{String,Any}(),
            "matrix_ctrl"   => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "matrix_ms"     => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "matrix_diff"   => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "sig_edges"     => Dict{String,Any}[],
            "band_stats"    => Dict{String,Any}[],
            "subject_means" => Dict{String,Any}[],
            "inclusion"     => Dict{String,Any}[],
        )

        gfe("transversal_summary.json") || return json(empty_resp)

        # ── Parse transversal_summary.json ────────────────────
        summ = Dict{String,Any}(
            "n_ms" => 0, "n_ctrl" => 0, "n_included" => 0,
            "n_excluded" => 0, "n_total_sig" => 0, "n_bands" => 0,
            "best_band" => "", "timestamp" => "",
        )
        try
            txt = read(joinpath(grp_dir, "transversal_summary.json"), String)
            for key in ["n_ms","n_ctrl","n_included","n_excluded","n_total_sig","n_bands"]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                m !== nothing && (summ[key] = parse(Int, m.captures[1]))
            end
            for key in ["best_band","cond","timestamp"]
                m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                m !== nothing && (summ[key] = m.captures[1])
            end
            # Band-level keys: ALPHA_n_sig, ALPHA_ctrl, etc.
            for m in eachmatch(r"\"([A-Z_]+)_(n_sig|ctrl|ms)\"\s*:\s*([0-9.]+)", txt)
                summ["$(m.captures[1])_$(m.captures[2])"] = tryparse(Float64, m.captures[3])
            end
        catch e
            @warn "transversal_summary.json parse error: $e"
        end

        # ── Función para leer matriz n×n ──────────────────────
        function read_matrix_csv(fname)
            p = joinpath(grp_dir, fname)
            isfile(p) || return Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}())
            try
                df  = CSV.read(p, DataFrame)
                ch  = string.(df[!, 1])
                n   = length(ch)
                vals = [[Float64(df[i, j+1]) for j in 1:n] for i in 1:n]
                Dict{String,Any}("channels" => ch, "values" => vals)
            catch
                Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}())
            end
        end

        matrix_ctrl = read_matrix_csv("group_connectivity_control_$(band).csv")
        matrix_ms   = read_matrix_csv("group_connectivity_ms_$(band).csv")
        matrix_diff = read_matrix_csv("group_difference_$(band).csv")

        # ── Significant edges ─────────────────────────────────
        sig_edges = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "significant_edges_$(band).csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(sig_edges, Dict{String,Any}(
                        "ch_a"     => string(row.ch_a),
                        "ch_b"     => string(row.ch_b),
                        "ctrl_mean"=> Float64(row.ctrl_mean),
                        "ms_mean"  => Float64(row.ms_mean),
                        "diff"     => Float64(row.diff),
                        "p_value"  => Float64(row.p_value),
                        "q_value"  => Float64(row.q_value),
                        "effect_d" => Float64(row.effect_d),
                    ))
                end
                sort!(sig_edges, by=r->abs(r["effect_d"]), rev=true)
            end
        catch e; @warn "significant_edges parse error: $e"; end

        # ── Band statistics ───────────────────────────────────
        band_stats = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "band_statistics.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(band_stats, Dict{String,Any}(
                        "band"       => string(row.band),
                        "n_channels" => Int(row.n_channels),
                        "n_pairs"    => Int(row.n_pairs),
                        "n_sig"      => Int(row.n_sig),
                        "pct_sig"    => Float64(row.pct_sig),
                        "ctrl_mean"  => Float64(row.ctrl_mean),
                        "ms_mean"    => Float64(row.ms_mean),
                        "diff_mean"  => Float64(row.diff_mean),
                        "mean_p"     => Float64(row.mean_p),
                        "mean_d"     => Float64(row.mean_d),
                    ))
                end
            end
        catch e; @warn "band_statistics parse error: $e"; end

        # ── Per-subject band means (for distribution) ─────────
        subject_means = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "subject_band_means.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(subject_means, Dict{String,Any}(
                        "subject_id" => string(row.subject_id),
                        "group"      => string(row.group),
                        "band"       => string(row.band),
                        "mean_wpli"  => Float64(row.mean_wpli),
                    ))
                end
            end
        catch e; @warn "subject_band_means parse error: $e"; end

        # ── Subject inclusion ─────────────────────────────────
        inclusion = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "subject_inclusion.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(inclusion, Dict{String,Any}(
                        "subject_id"      => string(row.subject_id),
                        "session_id"      => string(row.session_id),
                        "group"           => string(row.group),
                        "n_bands_ok"      => Int(row.n_bands_ok),
                        "included"        => Bool(row.included),
                        "excluded_reason" => string(row.excluded_reason),
                    ))
                end
            end
        catch e; @warn "subject_inclusion parse error: $e"; end

        json(Dict(
            "ok"            => true,
            "run"           => true,
            "cond"          => cond_short,
            "band"          => band,
            "summary"       => summ,
            "matrix_ctrl"   => matrix_ctrl,
            "matrix_ms"     => matrix_ms,
            "matrix_diff"   => matrix_diff,
            "sig_edges"     => sig_edges[1:min(50, end)],
            "band_stats"    => band_stats,
            "subject_means" => subject_means,
            "inclusion"     => inclusion,
        ))
    end

    # ─── API: Fase 14 — Evaluación Longitudinal ─────────────────
    route("/api/phase14_longitudinal") do
        cond_raw = string(get(getpayload(), :cond, "EC"))
        cond     = _normalize_cond(cond_raw)
        band     = uppercase(string(get(getpayload(), :band, "ALPHA")))
        cond_short = cond == "eyesclosed" ? "EC" : (cond == "eyesopen" ? "EO" : uppercase(cond_raw))
        grp_dir    = joinpath(res_root, "longitudinal", cond_short)

        gfe(f) = isfile(joinpath(grp_dir, f))

        empty_resp = Dict{String,Any}(
            "ok"            => true,
            "run"           => false,
            "cond"          => cond_short,
            "band"          => band,
            "summary"       => Dict{String,Any}(),
            "matrix_t1"     => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "matrix_t2"     => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "matrix_diff"   => Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()),
            "sig_edges"     => Dict{String,Any}[],
            "band_stats"    => Dict{String,Any}[],
            "subject_means" => Dict{String,Any}[],
            "paired_subjects"=> Dict{String,Any}[],
        )

        gfe("longitudinal_summary.json") || return json(empty_resp)

        summ = Dict{String,Any}(
            "n_paired"=>0,"n_t1"=>0,"n_t2"=>0,"n_loss"=>0,
            "n_total_sig"=>0,"n_bands"=>0,"best_band"=>"","timestamp"=>"",
        )
        try
            txt = read(joinpath(grp_dir, "longitudinal_summary.json"), String)
            for key in ["n_paired","n_t1","n_t2","n_loss","n_total_sig","n_bands"]
                m = match(Regex("\"$(key)\"\\s*:\\s*([0-9]+)"), txt)
                m !== nothing && (summ[key] = parse(Int, m.captures[1]))
            end
            for key in ["best_band","cond","timestamp"]
                m = match(Regex("\"$(key)\"\\s*:\\s*\"([^\"]+)\""), txt)
                m !== nothing && (summ[key] = m.captures[1])
            end
            for m in eachmatch(r"\"([A-Z_]+)_(n_sig|t1|t2)\"\s*:\s*([0-9.]+)", txt)
                summ["$(m.captures[1])_$(m.captures[2])"] = tryparse(Float64, m.captures[3])
            end
        catch e; @warn "longitudinal_summary.json parse error: $e"; end

        function read_matrix_csv(fname)
            p = joinpath(grp_dir, fname)
            isfile(p) || return Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}())
            try
                df  = CSV.read(p, DataFrame)
                ch  = string.(df[!, 1]); n = length(ch)
                vals = [[Float64(df[i, j+1]) for j in 1:n] for i in 1:n]
                Dict{String,Any}("channels"=>ch,"values"=>vals)
            catch; Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}()); end
        end

        matrix_t1   = read_matrix_csv("longitudinal_connectivity_t1_$(band).csv")
        matrix_t2   = read_matrix_csv("longitudinal_connectivity_t2_$(band).csv")
        matrix_diff = read_matrix_csv("longitudinal_difference_$(band).csv")

        sig_edges = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "significant_longitudinal_edges_$(band).csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(sig_edges, Dict{String,Any}(
                        "ch_a"    => string(row.ch_a),
                        "ch_b"    => string(row.ch_b),
                        "t1_mean" => Float64(row.t1_mean),
                        "t2_mean" => Float64(row.t2_mean),
                        "diff"    => Float64(row.diff),
                        "p_value" => Float64(row.p_value),
                        "q_value" => Float64(row.q_value),
                        "effect_d"=> Float64(row.effect_d),
                    ))
                end
                sort!(sig_edges, by=r->abs(r["effect_d"]), rev=true)
            end
        catch e; @warn "sig longitudinal edges: $e"; end

        band_stats = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "band_statistics_longitudinal.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(band_stats, Dict{String,Any}(
                        "band"      => string(row.band),
                        "n_channels"=> Int(row.n_channels),
                        "n_pairs"   => Int(row.n_pairs),
                        "n_sig"     => Int(row.n_sig),
                        "pct_sig"   => Float64(row.pct_sig),
                        "t1_mean"   => Float64(row.t1_mean),
                        "t2_mean"   => Float64(row.t2_mean),
                        "diff_mean" => Float64(row.diff_mean),
                        "mean_p"    => Float64(row.mean_p),
                        "mean_d"    => Float64(row.mean_d),
                    ))
                end
            end
        catch e; @warn "band stats longitudinal: $e"; end

        subject_means = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "subject_band_means.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(subject_means, Dict{String,Any}(
                        "subject_id"=> string(row.subject_id),
                        "timepoint" => string(row.timepoint),
                        "band"      => string(row.band),
                        "mean_wpli" => Float64(row.mean_wpli),
                    ))
                end
            end
        catch e; @warn "subject means longitudinal: $e"; end

        paired_subjects = Dict{String,Any}[]
        try
            p = joinpath(grp_dir, "paired_subjects.csv")
            if isfile(p)
                df = CSV.read(p, DataFrame)
                for row in eachrow(df)
                    push!(paired_subjects, Dict{String,Any}(
                        "subject_id"     => string(row.subject_id),
                        "session_t1"     => string(row.session_t1),
                        "session_t2"     => string(row.session_t2),
                        "n_bands_ok"     => Int(row.n_bands_ok),
                        "included"       => Bool(row.included),
                        "excluded_reason"=> string(row.excluded_reason),
                    ))
                end
            end
        catch e; @warn "paired subjects: $e"; end

        json(Dict(
            "ok"             => true,
            "run"            => true,
            "cond"           => cond_short,
            "band"           => band,
            "summary"        => summ,
            "matrix_t1"      => matrix_t1,
            "matrix_t2"      => matrix_t2,
            "matrix_diff"    => matrix_diff,
            "sig_edges"      => sig_edges[1:min(50, end)],
            "band_stats"     => band_stats,
            "subject_means"  => subject_means,
            "paired_subjects"=> paired_subjects,
        ))
    end

    # ─── Phase 15 — shared signal caches (used by validation + signal routes) ──
    _p15_json_cache = Dict{String, Any}()
    _p15_csv_cache  = Dict{String, Any}()

    function _p15_load_json!(path::String)
        haskey(_p15_json_cache, path) && return
        isfile(path) || return
        txt  = read(path, String)
        ch_m = match(r"\"channel_names\"\s*:\s*\[([^\]]+)\]", txt)
        chs  = ch_m !== nothing ?
            [strip(s, ['"',' ']) for s in split(ch_m.captures[1], ",")] : String[]
        tm_m = match(r"\"times\"\s*:\s*\[([^\]]+)\]", txt)
        times = tm_m !== nothing ?
            [parse(Float64, s) for s in split(tm_m.captures[1], ",")] : Float64[]
        entry = Dict{String,Any}("channel_names"=>chs, "times"=>times)
        for arr_key in ["raw_uv", "pre_ica_uv", "post_ica_uv", "filtered_uv"]
            kr = findfirst("\"$(arr_key)\"", txt)
            kr === nothing && continue
            outer = findnext('[', txt, last(kr))
            outer === nothing && continue
            pos = outer + 1
            ch_map = Dict{String, Vector{Float64}}()
            for ch in chs
                is = findnext('[', txt, pos)
                is === nothing && break
                depth = 0; ie = nothing
                for j in is:length(txt)
                    c = txt[j]
                    if c == '['; depth += 1
                    elseif c == ']'; depth -= 1; depth == 0 && (ie = j; break)
                    end
                end
                ie === nothing && break
                ch_map[ch] = [parse(Float64, strip(s))
                              for s in split(txt[is+1:ie-1], ",")]
                pos = ie + 1
            end
            entry[arr_key] = ch_map
        end
        _p15_json_cache[path] = entry
    end

    function _p15_load_csv!(path::String)
        haskey(_p15_csv_cache, path) && return
        isfile(path) || return
        try
            df    = CSV.read(path, DataFrame)
            times = Float64.(df[!, "t_s"])
            chs   = Dict{String, Vector{Float64}}()
            for col in names(df)
                col_str = string(col)
                col_str == "t_s" && continue
                chs[col_str] = Float64.(df[!, col])
            end
            _p15_csv_cache[path] = Dict{String,Any}("times"=>times, "channels"=>chs)
        catch e; @warn "p15 csv cache: $e" end
    end

    function _p15_sig_stats(v::Vector{Float64})::Dict{String,Float64}
        n = length(v); n == 0 && return Dict("rms"=>NaN,"std"=>NaN,"range"=>NaN)
        rms = sqrt(sum(v .^ 2) / n)
        μ   = sum(v) / n
        σ   = sqrt(max(0.0, sum((v .- μ) .^ 2) / (n - 1)))
        return Dict("rms"=>round(rms,digits=2), "std"=>round(σ,digits=2),
                    "range"=>round(maximum(v)-minimum(v),digits=2))
    end

    # ─── Phase 15 — MNE-Python cross-pipeline validation ─────
    route("/api/phase15_validation") do
        subj_raw   = string(get(getpayload(), :subj, ""))
        sess_raw   = string(get(getpayload(), :sess, ""))
        cond_raw   = string(get(getpayload(), :cond, "EC"))
        band       = uppercase(string(get(getpayload(), :band, "ALPHA")))
        cond_short = uppercase(cond_raw) in ["EC","EO"] ? uppercase(cond_raw) : "EC"
        task       = cond_short == "EC" ? "eyesclosed" : "eyesopen"
        mb_subj    = isempty(subj_raw) ? "M05" : subj_raw
        mb_sess    = isempty(sess_raw) ? "T2"  : sess_raw

        mb_root  = joinpath(project_root, "mne_brain")
        mb_res   = joinpath(mb_root, "results", "subjects",
                            "sub-$(mb_subj)", "ses-$(mb_sess)", task)
        val_dir  = joinpath(mb_root, "validation", "reports")
        tag      = "$(mb_subj)_$(mb_sess)"

        empty_resp = Dict{String,Any}(
            "ok"=>true, "run"=>false, "subj"=>mb_subj, "sess"=>mb_sess,
            "cond"=>cond_short, "band"=>band,
            "summary"=>Dict{String,Any}(), "matrix_julia"=>Dict{String,Any}(),
            "matrix_mne"=>Dict{String,Any}(), "figures"=>Dict{String,Any}())

        isdir(mb_res) || return json(empty_resp)

        # ── Parse comparison summary JSON via regex (no JSON3) ──────────
        summ     = Dict{String,Any}()
        summ_path = joinpath(val_dir, "comparison_$(tag)_$(cond_short)_summary.json")
        if isfile(summ_path)
            try
                txt = read(summ_path, String)
                for bk in ["DELTA","THETA","ALPHA","BETA_LOW","BETA_MID","BETA_HIGH","GAMMA"]
                    for mk in ["pearson_r","spearman_r","mae","rmse","bias","bias_pct","n"]
                        pat = Regex("\"$(bk)\"\\s*:\\s*\\{[^}]*\"$(mk)\"\\s*:\\s*([0-9.e+\\-]+)")
                        m   = match(pat, txt)
                        if m !== nothing
                            v = tryparse(Float64, m.captures[1])
                            v !== nothing && (summ["$(bk)_$(mk)"] = v)
                        end
                    end
                end
                for mk in ["psd_pearson_r","psd_mae_uv2","psd_bias_db_mean"]
                    m = match(Regex("\"$(mk)\"\\s*:\\s*([0-9.e+\\-]+)"), txt)
                    m !== nothing && (summ[mk] = tryparse(Float64, m.captures[1]))
                end
                for mk in ["nm_pipeline","mb_pipeline","note"]
                    m = match(Regex("\"$(mk)\"\\s*:\\s*\"([^\"]+)\""), txt)
                    m !== nothing && (summ[mk] = m.captures[1])
                end
            catch e; @warn "phase15 summary parse: $e"; end
        end

        # ── Read wPLI CSV → channels + 2D values matrix ─────────────────
        function _read_wpli(path::String)
            isfile(path) || return Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"vmax"=>0.0)
            try
                df   = CSV.read(path, DataFrame)
                chs  = string.(df[!, 1])
                n    = length(chs)
                vals = [[Float64(df[i, j+1]) for j in 1:n] for i in 1:n]
                vmax = isempty(vals) ? 0.0 : maximum(maximum.(vals))
                Dict{String,Any}("channels"=>chs,"values"=>vals,"vmax"=>vmax)
            catch e
                @warn "phase15 wpli read $path: $e"
                Dict("channels"=>String[],"values"=>Vector{Vector{Float64}}(),"vmax"=>0.0)
            end
        end

        julia_dir = joinpath(res_root, "subjects", "sub-$(mb_subj)", "ses-$(mb_sess)", task)
        mat_julia = _read_wpli(joinpath(julia_dir,  "wpli_$(band).csv"))
        mat_mne   = _read_wpli(joinpath(mb_res, "tables", "wpli_$(band)_$(cond_short).csv"))

        # ── Helper: CSV → array of row dicts ─────────────────────────────
        function _csv_rows(path::String)
            isfile(path) || return Vector{Dict{String,Any}}()
            try
                df = CSV.read(path, DataFrame)
                [Dict{String,Any}(string(c) => (ismissing(df[i,c]) ? nothing : df[i,c])
                                  for c in names(df))
                 for i in 1:nrow(df)]
            catch e; @warn "csv_rows $path: $e"; Vector{Dict{String,Any}}() end
        end

        # ── Helper: parse flat JSON keys (numerics + strings) ───────────
        function _json_keys(path::String, keys::Vector{String})
            isfile(path) || return Dict{String,Any}()
            try
                txt = read(path, String)
                out = Dict{String,Any}()
                for k in keys
                    m = match(Regex("\"$(k)\"\\s*:\\s*([0-9.e+\\-]+|true|false|\"[^\"]*\")"), txt)
                    m === nothing && continue
                    v = m.captures[1]
                    if v == "true";  out[k] = true
                    elseif v == "false"; out[k] = false
                    elseif startswith(v, "\""); out[k] = strip(v, '"')
                    else
                        n = tryparse(Float64, v)
                        n !== nothing && (out[k] = n)
                    end
                end
                out
            catch e; @warn "json_keys $path: $e"; Dict{String,Any}() end
        end

        # ── Helper: parse JSON array of objects (filter chain) ───────────
        function _json_array(path::String)
            isfile(path) || return Vector{Dict{String,Any}}()
            try
                txt  = read(path, String)
                objs = Vector{Dict{String,Any}}()
                for blk in eachmatch(r"\{([^}]+)\}", txt)
                    d = Dict{String,Any}()
                    for m2 in eachmatch(r"\"(\w+)\"\s*:\s*([0-9.e+\-]+|true|false|\"[^\"]*\")", blk.captures[1])
                        k2, v2 = m2.captures
                        if v2 == "true";  d[k2] = true
                        elseif v2 == "false"; d[k2] = false
                        elseif startswith(v2, "\""); d[k2] = strip(v2, '"')
                        else
                            n2 = tryparse(Float64, v2)
                            n2 !== nothing && (d[k2] = n2)
                        end
                    end
                    isempty(d) || push!(objs, d)
                end
                objs
            catch e; @warn "json_array $path: $e"; Vector{Dict{String,Any}}() end
        end

        # ── QC summaries ─────────────────────────────────────────────────
        qc_julia = _json_keys(
            joinpath(julia_dir, "artifact_rejection_summary.json"),
            ["n_total","n_valid","n_rejected","retention_pct","profile",
             "p2p_mean_uv","p2p_std_uv","p2p_max_uv","quality_mean",
             "quality_median","n_channels_used","n_rejected_amplitude",
             "n_rejected_gradient","min_amplitude_uv","max_amplitude_uv"])
        qc_mne = _json_keys(
            joinpath(mb_res, "tables", "epoch_summary_$(cond_short).json"),
            ["n_initial","n_valid","n_rejected","rejection_rate",
             "n_channels","epoch_duration_s","sfreq_hz"])

        # ── Filter chains ────────────────────────────────────────────────
        filter_mne = _json_array(joinpath(mb_res, "tables",
                                          "filter_chain_$(cond_short).json"))
        # NeuroMIND filter chain is fixed (from pipeline.toml)
        filter_julia = [
            Dict{String,Any}("step"=>1,"name"=>"Notch",
                "freq"=>"49.5–50.5 Hz","method"=>"filt","order"=>4,"applied"=>true),
            Dict{String,Any}("step"=>2,"name"=>"Bandreject",
                "freq"=>"99.5–100.5 Hz","method"=>"filt","order"=>4,"applied"=>true),
            Dict{String,Any}("step"=>3,"name"=>"High-pass",
                "freq"=>"0.5 Hz","method"=>"filtfilt","order"=>4,"applied"=>true),
            Dict{String,Any}("step"=>4,"name"=>"Low-pass",
                "freq"=>"150.0 Hz","method"=>"filtfilt","order"=>4,"applied"=>true),
        ]

        # ── Per-channel stats ────────────────────────────────────────────
        ch_stats_julia = _csv_rows(joinpath(julia_dir, "channel_statistics.csv"))
        ch_stats_mne   = _csv_rows(joinpath(mb_res, "tables",
                                            "qc_channels_$(cond_short).csv"))

        # channel_artifact_summary → bad epochs per channel (Julia only)
        ch_art_rows = _csv_rows(joinpath(julia_dir, "channel_artifact_summary.csv"))
        ch_art_julia = Dict{String,Any}()
        for r in ch_art_rows
            ch = string(get(r, "channel", ""))
            isempty(ch) || (ch_art_julia[ch] = r)
        end

        # ── Per-channel band power ────────────────────────────────────────
        function _read_band_power(path::String)
            isfile(path) || return Dict{String,Any}()
            try
                df  = CSV.read(path, DataFrame)
                chs = string.(df[!, 1])
                cols = [string(c) for c in names(df)[2:end]]
                Dict{String,Any}(chs[i] =>
                    Dict{String,Any}(cols[j] => Float64(df[i, j+1])
                                     for j in eachindex(cols)
                                     if !ismissing(df[i, j+1]))
                    for i in eachindex(chs))
            catch e; @warn "band_power $path: $e"; Dict{String,Any}() end
        end
        bp_julia = _read_band_power(joinpath(julia_dir, "band_power_summary.csv"))
        bp_mne   = _read_band_power(joinpath(mb_res, "tables",
                                             "band_power_summary_$(cond_short).csv"))

        # ── 6-way channel statistics (raw + pre-ICA + post-ICA, Julia + MNE) ───
        sig_preview_path  = joinpath(mb_res, "tables", "signal_preview_$(cond_short).json")
        julia_raw_path    = joinpath(julia_dir, "raw_signal.csv")
        julia_preica_path = joinpath(julia_dir, "ica_signal_before.csv")
        julia_ica_path    = joinpath(julia_dir, "ica_signal_after.csv")
        _p15_load_json!(sig_preview_path)
        _p15_load_csv!(julia_raw_path)
        _p15_load_csv!(julia_preica_path)
        _p15_load_csv!(julia_ica_path)

        ch_stats_6way = Dict{String, Any}()
        # Julia 3 stages from CSV files
        for (lbl, cpath) in [("julia_raw", julia_raw_path),
                              ("julia_preica", julia_preica_path),
                              ("julia_ica", julia_ica_path)]
            haskey(_p15_csv_cache, cpath) || continue
            for (ch, v) in _p15_csv_cache[cpath]["channels"]
                haskey(ch_stats_6way, ch) || (ch_stats_6way[ch] = Dict{String,Any}())
                ch_stats_6way[ch][lbl] = _p15_sig_stats(v)
            end
        end
        # MNE 3 stages from preview JSON
        if haskey(_p15_json_cache, sig_preview_path)
            jc = _p15_json_cache[sig_preview_path]
            for (arr_key, lbl) in [("raw_uv","mne_raw"),
                                    ("pre_ica_uv","mne_preica"),
                                    ("post_ica_uv","mne_ica")]
                arr = get(jc, arr_key, Dict{String,Vector{Float64}}())
                # backward-compat: old JSONs used filtered_uv for pre-ICA
                if isempty(arr) && lbl in ("mne_preica","mne_ica")
                    arr = get(jc, "filtered_uv", Dict{String,Vector{Float64}}())
                end
                for (ch, v) in arr
                    haskey(ch_stats_6way, ch) || (ch_stats_6way[ch] = Dict{String,Any}())
                    ch_stats_6way[ch][lbl] = _p15_sig_stats(v)
                end
            end
        end
        # keep old key for backward compat with any cached dashboard state
        ch_stats_4way = ch_stats_6way

        # ── Base64-encode PNG figures ────────────────────────────────────
        figs = Dict{String,String}()
        for (k, fname) in [
            ("wpli",       "comparison_$(tag)_$(cond_short)_wpli.png"),
            ("heatmaps",   "comparison_$(tag)_$(cond_short)_heatmaps.png"),
            ("psd",        "comparison_$(tag)_$(cond_short)_psd.png"),
            ("band_power", "comparison_$(tag)_$(cond_short)_band_power.png"),
        ]
            p = joinpath(val_dir, fname)
            isfile(p) && (figs[k] = base64encode(read(p)))
        end

        return json(Dict{String,Any}(
            "ok"             => true,
            "run"            => !isempty(mat_mne["channels"]),
            "subj"           => mb_subj,
            "sess"           => mb_sess,
            "cond"           => cond_short,
            "band"           => band,
            "summary"        => summ,
            "matrix_julia"   => mat_julia,
            "matrix_mne"     => mat_mne,
            "figures"        => figs,
            "qc_julia"       => qc_julia,
            "qc_mne"         => qc_mne,
            "filter_julia"   => filter_julia,
            "filter_mne"     => filter_mne,
            "ch_stats_julia"  => ch_stats_julia,
            "ch_stats_mne"    => ch_stats_mne,
            "ch_art_julia"    => ch_art_julia,
            "ch_stats_4way"   => ch_stats_4way,
            "bp_julia"        => bp_julia,
            "bp_mne"          => bp_mne,
        ))
    end

    # ─── Phase 15: signal comparison endpoint ─────────────────
    route("/api/phase15_signal") do
        subj_raw = string(get(getpayload(), :subj, "M05"))
        sess_raw = string(get(getpayload(), :sess, "T2"))
        cond_raw = string(get(getpayload(), :cond, "EC"))
        channel  = string(get(getpayload(), :channel, "Fz"))
        mode_raw = lowercase(string(get(getpayload(), :mode, "raw")))
        # 3 valid modes: raw | pre_ica | post_ica
        mode = mode_raw in ["pre_ica","post_ica"] ? mode_raw : "raw"
        cond_s   = uppercase(cond_raw) in ["EC","EO"] ? uppercase(cond_raw) : "EC"
        task     = cond_s == "EC" ? "eyesclosed" : "eyesopen"
        mb_root  = joinpath(project_root, "mne_brain")
        julia_dir = joinpath(res_root, "subjects", "sub-$(subj_raw)", "ses-$(sess_raw)", task)
        mb_res   = joinpath(mb_root, "results", "subjects",
                            "sub-$(subj_raw)", "ses-$(sess_raw)", task)

        # ── NeuroMIND (Julia) ─────────────────────────────────
        julia_sig = Dict{String,Any}()
        sig_file, julia_lbl = if mode == "raw"
            "raw_signal.csv", "NeuroMIND (Julia) — cruda (sin filtrar)"
        elseif mode == "pre_ica"
            "ica_signal_before.csv", "NeuroMIND (Julia) — filtrada pre-ICA"
        else
            "ica_signal_after.csv", "NeuroMIND (Julia) — post-ICA"
        end
        sig_path  = joinpath(julia_dir, sig_file)
        _p15_load_csv!(sig_path)
        if haskey(_p15_csv_cache, sig_path)
            c = _p15_csv_cache[sig_path]
            if haskey(c["channels"], channel)
                julia_sig = Dict{String,Any}(
                    "times"=>c["times"], "signal"=>c["channels"][channel],
                    "fs"=>500.0, "n_samples"=>length(c["channels"][channel]),
                    "label"=>julia_lbl)
            end
        end

        # ── mne_brain (MNE-Python) ────────────────────────────
        mne_sig   = Dict{String,Any}()
        prev_path = joinpath(mb_res, "tables", "signal_preview_$(cond_s).json")
        mne_arr, mne_lbl = if mode == "raw"
            "raw_uv", "mne_brain (MNE) — cruda (sin filtrar)"
        elseif mode == "pre_ica"
            # prefer new key, fall back to old filtered_uv
            "pre_ica_uv", "mne_brain (MNE) — filtrada pre-ICA"
        else
            # prefer new key, fall back to old filtered_uv
            "post_ica_uv", "mne_brain (MNE) — post-ICA"
        end
        _p15_load_json!(prev_path)
        if haskey(_p15_json_cache, prev_path)
            c = _p15_json_cache[prev_path]
            arr = get(c, mne_arr, Dict{String,Vector{Float64}}())
            # fallback for old JSONs that used filtered_uv for both pre- and post-ICA
            if isempty(arr) && mode == "pre_ica"
                arr = get(c, "filtered_uv", Dict{String,Vector{Float64}}())
            elseif isempty(arr) && mode == "post_ica"
                arr = get(c, "filtered_uv", Dict{String,Vector{Float64}}())
            end
            if haskey(arr, channel)
                mne_sig = Dict{String,Any}(
                    "times"=>c["times"], "signal"=>arr[channel],
                    "fs"=>500.0, "n_samples"=>length(arr[channel]),
                    "label"=>mne_lbl)
            end
        end

        return json(Dict{String,Any}(
            "ok"=>true, "channel"=>channel, "cond"=>cond_s, "mode"=>mode,
            "julia"=>julia_sig, "mne"=>mne_sig))
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
