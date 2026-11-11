# MINDEM-IMIBIC Resting-State EEG

Resting-state EEG of patients with multiple sclerosis (MS) and healthy
controls (HC), recorded in two conditions (`task-eyesclosed`,
`task-eyesopen`) and, for MS patients, at two time points
(`ses-T1` baseline, `ses-T2` follow-up after an intervention of about
12 weeks; 77-119 days between assessments).
Healthy controls have a single session (`ses-T1`). Controls MC5, MC6, MC11
and MC12 were excluded by the clinical team because their EEG was too
noisy; MC39 and MC40 were recruited later (the control numbering is not
contiguous).

This dataset is generated, not edited by hand: `scripts/audit_full_dataset.jl`
writes `participants.tsv` and `scripts/build_bids_full.jl` writes everything
else from `data/sourcedata/` (original BrainVision recordings and clinical
spreadsheets) following the rules in `config/raw_selection.csv`.

## Acquisition

- 31 EEG channels, international 10-20 system (ECI cap), 500 Hz,
  BrainVision Recorder 1.21 with a Brain Products actiCHamp amplifier.
- Online reference FCz. Ground: dedicated frontal `Gnd` electrode of the cap.
- Analogue filters (amplifier): DC to 140 Hz. Recorder software filters:
  time constant 1.59155 s (about 0.1 Hz) high-pass, 70 Hz low-pass, 50 Hz notch.
  They are already applied to the stored signal.
- Signals are stored as BrainVision `IEEE_FLOAT_32`, `MULTIPLEXED`; the
  channel resolution (0.0488281 uV) in the `.vhdr` must be applied to obtain
  microvolts.

## Recordings cut from a single raw file

For sub-M04, sub-M05 and sub-M06 at `ses-T1` the eyes-open and eyes-closed
blocks were recorded in one file. They were split before any processing at
the seconds chosen by the clinical team, with half-open intervals
`[tmin, tmax)`, so that the sample at the transition belongs to one block only.
Each of the six resulting recordings has its own sidecar
`*_eeg.json` with the `SourceWindow` used. The signal is otherwise
unmodified.

## Participants

`participants.tsv` includes clinical variables (EDSS, D-FIS, BDI-II,
diagnosis, treatment) and dates of birth and inclusion. Column definitions
are in `participants.json`. Fatigue categories (`dfis_group_*`,
`fatigue_response`) split each sample at its own D-FIS median, with ties
going to the upper category; the clinical team's original categories are
kept in the `*_orig` columns. Primary analyses use continuous D-FIS. **The dataset contains personal health data and
must not be shared or published without removing the identifying columns.**

## Validation

`bids-validator-deno data/bids` (validator 3.0.2) reports no errors. The
remaining warnings are known and accepted:

- `EVENTS_TSV_MISSING`: resting-state task without stimulus events; the
  BrainVision markers (`*.vmrk`) only contain `New Segment` entries.
- `TOO_FEW_AUTHORS`: single author listed.
- `JSON_KEY_RECOMMENDED` (`HEDVersion`, `SourceDatasets`): no HED tags are
  used, and the source data (`data/sourcedata/`) is not a BIDS dataset.
- `SIDECAR_KEY_RECOMMENDED` (`Instructions`): the task instructions given to
  the participants are not documented in the source protocol files.

## Electrode coordinates

`*_electrodes.tsv` come from a normalised cap template; their origin and axes
are not documented in the source, so the coordinate system is declared as
`Other` in `*_coordsystem.json`.
