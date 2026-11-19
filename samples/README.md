# Audio samples

English · [日本語](README.ja.md)

Each WAV is 48 kHz, 16-bit, mono. The full performances are ascending B-flat major scales at concert pitch with A4 = 440 Hz.

## Three-octave scales

| Sample | Concert range | Duration | Regenerate from the repository root |
| --- | --- | --- | --- |
| [York trumpet with Bach 5C](trumpet_Bach5C_Bb3_Bb6_major_f_strong_attack.wav) | B♭3–B♭6 | 26.31 s | `python src/render_trumpet_bach5c.py` |
| [Conn 8H trombone with Wick 5AL](trombone_5AL_Bb2_Bb5_major_f_strong_attack.wav) | B♭2–B♭5 | 26.31 s | `python src/render_trombone_detailed_scale.py` |
| [5AL trombone, louder forte](trombone_5AL_Bb2_Bb5_major_louder_f_strong_attack.wav) | B♭2–B♭5 | 26.31 s | `python src/render_trombone_detailed_scale.py --louder` |
| [5AL trombone, stronger attack](trombone_5AL_Bb2_Bb5_major_stronger_f_strong_attack.wav) | B♭2–B♭5 | 26.31 s | `python src/render_trombone_detailed_scale.py --stronger` |
| [BB♭ tuba with Yamaha BB-67C4](BBb_tuba_Bb1_Bb4_major_f_strong_attack.wav) | B♭1–B♭4 | 26.31 s | `python src/render_tuba_scale.py` |

The three-octave scales have 22 notes. Each note starts from rest and runs for 1.08 s, with 1.65 s for the final note and 90 ms of silence after each note. Generate the three trombone variants in the order shown so the louder and stronger runs can use the preceding variant's output gain and calculate the acoustic level changes.

## Trombone listening copies

| Sample | Duration | Processing |
| --- | --- | --- |
| [Stronger attack, peak-limited listening copy](trombone_5AL_Bb2_Bb5_major_stronger_f_listening.wav) | 26.31 s | Playback copy with peak limiting. Level and peak measurements are in the [stronger-attack report](../docs/results/trombone_5AL_Bb2_Bb5_stronger_f_measurements.json). |
| [Louder and stronger attacks, matched gain](trombone_5AL_stronger_forte_before_after.wav) | 9.18 s | Each pair plays the louder preset followed by the stronger preset. The pairs are B♭2, B♭4, and B♭5 in that order. Their pressure scales are restored before one gain is applied to the assembled WAV. |

The full simulation WAVs use a common gain within each file. For acoustic comparisons between files, account for their different gains using the [measurement reports](../docs/results/).

## Earlier scales

| File | Duration | Regenerate from the repository root |
| --- | --- | --- |
| [Trumpet high-register reference](trumpet_Bb3_Bb6_major_f_strong_attack.wav) | 26.31 s | `python src/render_trumpet_high.py` |
| [Two-octave trumpet with thermoviscous relaxation](trumpet_Bb3_Bb5_major_f_modern_1p5x.wav) | 18.12 s | `python src/render_trumpet_modern.py` |
| [Trumpet with linear bore](trumpet_Bb3_Bb5_major_f.wav) | 12.08 s | `python src/render_scale.py --render` |
| [Trumpet with prescribed-wave propagation](trumpet_Bb3_Bb5_major_f_nonlinear.wav) | 12.08 s | `python src/render_nonlinear.py` |
| [Trumpet with bidirectional waveguide](trumpet_Bb3_Bb5_major_f_bidirectional_nonlinear.wav) | 12.08 s | `python src/render_bidirectional.py` |
| [Original axisymmetric trumpet](trumpet_Bb3_Bb5_major_f_axisymmetric.wav) | 12.08 s | `python src/render_axisymmetric.py` |
| [Standard trombone](trombone_Bb2_Bb4_major_f_axisymmetric.wav) | 12.08 s | `python src/trombone_axisymmetric.py` |
| [Trombone with stronger attack](trombone_Bb2_Bb4_major_f_strong_attack.wav) | 12.08 s | `python src/render_trombone_attack.py --preset data/performance_trombone_attack.json` |

Generate the earlier standard trombone performance before its two-octave stronger-attack version so the comparison can use the standard measurements and traces.

The original 12.08 s scales use 0.72 s notes, a 1.10 s final note, and 0.06 s gaps. The 18.12 s trumpet version uses 1.08 s notes, a 1.65 s final note, and 0.09 s gaps.

[The manifest](manifest.json) records each sample's format, duration, and SHA-256 checksum. Simulations on other drivers or hardware may differ slightly in floating-point details.
