# Audio samples

Each WAV is 48 kHz, 16-bit, mono. All performances use an ascending B-flat major scale at concert pitch with A4 = 440 Hz. The trumpet covers B♭3–B♭5 and the tenor trombone covers B♭2–B♭4.

| File | Duration | Regenerate from the repository root |
| --- | --- | --- |
| [Current trumpet](trumpet_Bb3_Bb5_major_f_modern_1p5x.wav) | 18.12 s | `python src/render_trumpet_modern.py` |
| [Trumpet with linear bore](trumpet_Bb3_Bb5_major_f.wav) | 12.08 s | `python src/render_scale.py --render` |
| [Trumpet with prescribed-wave propagation](trumpet_Bb3_Bb5_major_f_nonlinear.wav) | 12.08 s | `python src/render_nonlinear.py` |
| [Trumpet with bidirectional waveguide](trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear.wav) | 12.08 s | `python src/render_bidirectional.py` |
| [Original axisymmetric trumpet](trumpet_Bb3_Bb5_major_f_axisymmetric.wav) | 12.08 s | `python src/render_axisymmetric.py` |
| [Standard trombone](trombone_Bb2_Bb4_major_f_axisymmetric.wav) | 12.08 s | `python src/trombone_axisymmetric.py` |
| [Trombone with stronger attack](trombone_Bb2_Bb4_major_f_strong_attack.wav) | 12.08 s | `python src/render_trombone_attack.py --preset data/performance_trombone_attack.json` |

Generate the standard trombone performance before the stronger-attack version, which uses the standard measurements and traces for its comparison.

The original scales use 0.72 s notes, a 1.10 s final note, and 0.06 s gaps, including a final gap. The current trumpet uses 1.08 s notes, a 1.65 s final note, and 0.09 s gaps. Each note is simulated from rest.

The audio is computed from the modeled radiated pressure. Each complete file receives a common gain before conversion to PCM. File gains differ between model versions. A level-matched comparison requires accounting for the recorded gains in the corresponding [measurement reports](../docs/results/).

The files were generated on 17 November 2026. `manifest.json` records each sample's format, duration, and SHA-256 checksum. Simulations on other drivers or hardware may differ slightly in floating-point details.
