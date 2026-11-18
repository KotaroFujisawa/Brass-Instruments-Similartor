# Geometry, presets, and provenance

All lengths are in metres and all frequencies are in hertz unless a field name states another unit. Bore tables describe distance along an unfolded instrument centerline. The source graphs and the inferred dimensions are documented in [the model notes](../docs/models.md) and [references](../docs/references.md).

| File | Contents |
| --- | --- |
| `york_trumpet_fig7_52.csv` | Distance and diameter digitized from the York trumpet in Figure 7.52 of *The Science of Brass Instruments* |
| `conn8h_fig4_22.csv` | Distance and radius of the reconstructed Conn 8H body, beginning 80 mm from the mouthpiece entrance |
| `wick5al_calibrated.csv` | Distance and radius of the calibrated representative 5AL mouthpiece |
| `model.json` | Trumpet cup depth and the original detached-mouthpiece calibration target |
| `performance.json` | Linear-bore trumpet performance |
| `performance_bidirectional.json` | Bidirectional CPU trumpet performance |
| `performance_axisymmetric.json` | Original axisymmetric trumpet performance |
| `performance_trumpet_modern.json` | Per-note overrides applied to the original axisymmetric trumpet preset by the current renderer |
| `performance_trombone.json` | Trombone notes, slide positions, lip settings, blowing pressures, and length corrections |
| `performance_trombone_attack.json` | Trombone performance with selected attack-envelope settings |
| `sources.json` | Provenance of the initial trumpet geometry, mouthpiece reconstruction, and lip model |
| `nonlinear_sources.json` | Equations and assumptions for prescribed-wave propagation |
| `sources_bidirectional.json` | References and assumptions for the bidirectional CPU model |

Performance fields include `target_hz`, `lip_frequency_hz`, `blowing_pressure_pa`, and `rest_y_m`. The trumpet's `valves_semitones` lists the nominal lowering of each active valve. `tuning_m` and the trombone's `trim_m` adjust effective cylindrical bore length for the numerical performance.

The bore profiles are approximate graph reconstructions. Mouthpiece lengths, volumes, and resonances constrain inferred shapes. Performance presets are fitted playing settings, and their values should be interpreted together with the model version.
