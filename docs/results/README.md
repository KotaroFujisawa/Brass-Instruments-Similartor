# Published measurements

These reports accompany the committed audio samples. They were generated during the model-development runs on 17 November 2026. Re-running the renderers and verification scripts writes new reports under `output/`.

| Report | Contents |
| --- | --- |
| [Current trumpet measurements](trumpet_modern_measurements.json) | Per-note pitches, early and late pitch estimates, pressure, periodicity, and solver settings |
| [Current trumpet verification](trumpet_modern_verification.json) | Weak-pulse impedance, bell reflection magnitude, and a B♭4 feedback diagnostic |
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
