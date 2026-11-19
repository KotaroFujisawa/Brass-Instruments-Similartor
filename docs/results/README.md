# Published measurements

English · [日本語](README.ja.md)

These reports accompany the [published audio samples](../../samples/). Re-running the renderers and verification scripts writes new reports under `output/`. Sample paths in the published reports are relative to the repository root.

## Three-octave performances and model checks

| Report | Contents |
| --- | --- |
| [Trumpet with Bach 5C](trumpet_Bach5C_Bb3_Bb6_measurements.json) | Per-note render data and independent analysis of the 22-note PCM file |
| [Trumpet high-register reference](trumpet_high_measurements.json) | Earlier B♭3–B♭6 reference used in the Bach 5C comparison |
| [Trombone with Wick 5AL](trombone_5AL_Bb2_Bb5_measurements.json) | Per-note render data and independent B♭2–B♭5 WAV analysis |
| [Louder 5AL trombone](trombone_5AL_Bb2_Bb5_louder_f_measurements.json) | Paired acoustic and PCM level measurements relative to the detailed 5AL scale |
| [5AL trombone with stronger attack](trombone_5AL_Bb2_Bb5_stronger_f_measurements.json) | Onset-aligned and fixed-time attack changes, sustain levels, and peak-limited listening-copy measurements |
| [BB♭ tuba with Yamaha BB-67C4](BBb_tuba_Bb1_Bb4_measurements.json) | B♭1–B♭4 render and independent PCM analysis, including attack timing and early low-note modulation |
| [Trumpet taper and bell comparison](bell_comparison.json) | Comparisons of the earlier and refined tapers and the coupled exterior radiation model. [Taper plot](../figures/trumpet_taper_refinement.png) |
| [Bell boundary check](bell_boundary_verification.json) | Comparison of CUDA boundary pressures with the CPU frequency-domain reconstruction |
| [Bell exterior verification](bell_exterior_verification.json) | Aperture projection, state-space fit, passivity, and receiver transfer checks |
| [Trumpet mouthpiece comparison](trumpet_mouthpiece_results.json) | Geometry, impedance, and playing tests comparing the earlier mouthpiece and the measured Bach 5C. [Impedance plot](../figures/mouthpiece_impedance.png) |

## Earlier results

| Report | Contents |
| --- | --- |
| [Two-octave trumpet measurements](trumpet_modern_measurements.json) | Per-note pitches, early and late pitch estimates, pressure, periodicity, and solver settings |
| [Two-octave trumpet verification](trumpet_modern_verification.json) | Weak-pulse impedance, bell reflection magnitude, and a B♭4 feedback diagnostic |
| [Euler solver tests](axisymmetric_solver_tests.json) | Rest, conservation, Sod shock-tube refinement, and radial-mode checks |
| [Original axisymmetric trumpet](trumpet_axisymmetric_measurements.json) | Per-note results for the earlier axisymmetric sample |
| [Trombone measurements](trombone_measurements.json) | Standard performance, presets, and per-note diagnostics |
| [Trombone geometry](trombone_geometry_verification.json) | Mouthpiece volume and resonance, slide lengths, and impedance peaks |
| [Trombone coupling](trombone_coupling_tests.json) | B♭3 feedback diagnostic and the tested mesh refinement |
| [Stronger trombone attack](trombone_strong_attack_measurements.json) | Attack settings, onset comparisons, pitch, and level changes |
| [Linear-bore trumpet](trumpet_linear_measurements.json) | Original CPU performance |
| [Bidirectional CPU trumpet](trumpet_bidirectional_measurements.json) | Nonlinear and linear comparison metrics |
| [Prescribed-wave trumpet](trumpet_one_way_measurements.json) | Nonlinear propagation and matched-reference metrics |

Stored computation times depend on hardware and concurrent work during the development runs. The `common_gain_fs_per_pa`, `gain_fs_per_pa`, and related gain fields describe PCM scaling for their respective renderers.
