# Geometry, presets, and provenance

All lengths are in metres and all frequencies are in hertz unless a field name states another unit. Bore tables describe distance along an unfolded instrument centerline. The source graphs and the inferred dimensions are documented in [the model notes](../docs/models.md) and [references](../docs/references.md).

| File | Contents |
| --- | --- |
| `york_trumpet_fig7_52.csv` | Distance and diameter digitized from the York trumpet in Figure 7.52 of *The Science of Brass Instruments* |
| `bach_5c_jackson_2025_fig2.csv` | Axial distance and radius digitized from the measured Bach 5C contour in Jackson and Scavone (2025), Figure 2 |
| `conn8h_fig4_22.csv` | Distance and radius of the reconstructed Conn 8H body, beginning 80 mm from the mouthpiece entrance |
| `wick5al_calibrated.csv` | Distance and radius of the calibrated representative 5AL mouthpiece |
| `alexander_bb_tuba_reconstructed.csv` | Tuba bore reconstructed for BB♭ from the Alexander C-tuba profile, before the solver's 0.936 axial calibration |
| `yamaha_bb67c4_reconstructed.csv` | BB-67C4 internal profile reconstructed from measured dimensions and volumes |
| `model.json` | Trumpet cup depth and the original detached-mouthpiece calibration target |
| `performance.json` | Linear-bore trumpet performance |
| `performance_bidirectional.json` | Bidirectional CPU trumpet performance |
| `performance_axisymmetric.json` | Original axisymmetric trumpet performance |
| `performance_trumpet_modern.json` | Per-note overrides for the earlier two-octave trumpet with cubic stiffness and thermoviscous losses |
| `performance_trumpet_high.json` | Extended B♭3–B♭6 trumpet preset and attack settings |
| `performance_trumpet_bach5c.json` | Trumpet overrides for the refined Bach 5C geometry and coupled exterior |
| `performance_trombone.json` | Trombone notes, slide positions, lip settings, blowing pressures, and length corrections |
| `performance_trombone_attack.json` | Trombone performance with selected attack-envelope settings |
| `performance_trombone_detailed_scale.json` | Overrides and high-register settings for the three-octave 5AL trombone |
| `performance_trombone_stronger_forte.json` | Tuning overrides for the `--stronger` trombone performance |
| `performance_tuba.json` | Tuba lip settings, pressures, attack envelopes, valve selections, and length offsets for B♭1–B♭4 |
| `sources.json` | Provenance of the initial trumpet geometry, mouthpiece reconstruction, and lip model |
| `nonlinear_sources.json` | Equations and assumptions for prescribed-wave propagation |
| `sources_bidirectional.json` | References and assumptions for the bidirectional CPU model |

Performance fields include `target_hz`, `lip_frequency_hz`, `blowing_pressure_pa`, and `rest_y_m`. The trumpet's `valves_semitones` lists the nominal lowering of each active valve. The tuba's `valves` lists valve numbers 1–4, whose nominal lowerings are two, one, three, and five semitones in that order. `tuning_m` and `trim_m` adjust effective bore length. The trombone's `--louder` pressure multiplier is defined in the renderer.

The bore profiles are approximate graph reconstructions. Measured contours and available mouthpiece lengths, volumes, and resonances constrain the selected geometry. Performance presets are fitted to their model version.

`bach_5c_jackson_2025_fig2.csv` reproduces the measured contour from Figure 2 by sampling and converting its published coordinate path. Credit Miranda Jackson and Gary Scavone, “A study of impedance of brass instruments and mouthpieces,” *Proceedings of Meetings on Acoustics* 58, 035019 (2025). The source article is published under CC BY-NC 4.0. [DOI and source figure](https://doi.org/10.1121/2.0002165) · [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/)
