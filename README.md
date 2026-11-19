# Brass Instrument Simulator

Physical sound synthesis for B-flat trumpet, tenor trombone, and BB-flat contrabass tuba. The models couple the two-coordinate lip motion of Adachi and Sato to a mouthpiece and a bore. The CUDA solver models nonlinear propagation in the unfolded axisymmetric bore and couples the bell to a reduced exterior sound field. Reflections from the bell travel back through the bore and act on the lips.

[Model equations and assumptions](docs/models.md) · [References](docs/references.md) · [Geometry and presets](data/README.md) · [Japanese notes](docs/ja/README.md)

## Audio samples

These forte scales run from low to high in B-flat major at concert pitch (A4 = 440 Hz). Each contains 22 independently simulated notes and lasts 26.31 s. WAV files are 48 kHz, 16-bit mono. Most notes last 1.08 s and the final high B-flat lasts 1.65 s, with a 90 ms gap after each note.

| Sample | Range | Measurements |
| --- | --- | --- |
| [York trumpet with Bach 5C](samples/trumpet_Bach5C_Bb3_Bb6_major_f_strong_attack.wav) | B♭3–B♭6 | [Pitch and radiated pressure](docs/results/trumpet_Bach5C_Bb3_Bb6_measurements.json) |
| [Conn 8H trombone with Denis Wick 5AL](samples/trombone_5AL_Bb2_Bb5_major_f_strong_attack.wav) | B♭2–B♭5 | [Base forte](docs/results/trombone_5AL_Bb2_Bb5_measurements.json) |
| [Trombone with higher blowing pressure](samples/trombone_5AL_Bb2_Bb5_major_louder_f_strong_attack.wav) | B♭2–B♭5 | [15% increase in steady blowing pressure](docs/results/trombone_5AL_Bb2_Bb5_louder_f_measurements.json) |
| [Trombone with higher pressure and a stronger attack](samples/trombone_5AL_Bb2_Bb5_major_stronger_f_strong_attack.wav) | B♭2–B♭5 | [Pressure and onset comparison](docs/results/trombone_5AL_Bb2_Bb5_stronger_f_measurements.json) |
| [Reconstructed BB♭ tuba with Yamaha BB-67C4](samples/BBb_tuba_Bb1_Bb4_major_f_strong_attack.wav) | B♭1–B♭4 | [Pitch, pressure, and onset](docs/results/BBb_tuba_Bb1_Bb4_measurements.json) |

Relative to the higher-pressure trombone sample, the stronger-attack version adds another 25% to steady blowing pressure and reduces its rise time from 12 to 6.5 ms. The simulated median pressure level rises by 3.8 dB during sustain and by 6.0 dB in the first 50 ms after acoustic onset. The main scale files use an individual gain for each complete file. A [common-gain excerpt](samples/trombone_5AL_stronger_forte_before_after.wav) preserves the level change, and a [peak-limited listening copy](samples/trombone_5AL_Bb2_Bb5_major_stronger_f_listening.wav) keeps the strongest attacks from setting playback level for the entire scale. [Earlier models, checksums, and sample details](samples/README.md)

## Install

Python 3.12 was used for development. NumPy, SciPy, Numba, and Matplotlib are required. Run commands from the repository root.

```sh
git clone https://github.com/KotaroFujisawa/Brass-Instruments-Similartor.git
cd Brass-Instruments-Similartor
python -m venv .venv
```

Activate the environment on Windows PowerShell with `.venv\Scripts\Activate.ps1`, or on Linux with `source .venv/bin/activate`, then install the dependencies.

```sh
python -m pip install -r requirements.txt
```

The CPU models can run with these dependencies. The axisymmetric models also require an NVIDIA GPU, an NVIDIA driver, and NVRTC, the CUDA runtime compiler. The solver uses the CUDA Driver API directly through Python ctypes.

```sh
python -m pip install nvidia-cuda-nvrtc-cu12==12.4.127
```

The CUDA calculations were run on an NVIDIA RTX 2000 Ada Generation with NVRTC 12.4.127. On Windows, the loader searches installed Python packages, `CUDA_NVRTC_DIR`, and `CUDA_PATH/bin`. On Linux, make `libnvrtc.so.12` and its companion libraries available through the system loader or `LD_LIBRARY_PATH`. Linux execution has not been verified in this project.

## Generate the three-octave scales

```sh
python src/render_trumpet_bach5c.py
python src/render_trombone_detailed_scale.py
python src/render_trombone_detailed_scale.py --louder
python src/render_trombone_detailed_scale.py --stronger
python src/render_tuba_scale.py
```

The WAV and per-note measurements are written under `output/trumpet_bach5c/`, `output/trombone_detailed_scale/`, `output/trombone_louder_forte/`, `output/trombone_stronger_forte/`, and `output/tuba_axisymmetric/`, respectively. Generate the trombone variants in the listed order to include their pressure comparisons to the preceding version. The presets specify the lip frequency, blowing pressure, rest opening, fingerings or slide positions, and length offsets for each note.

To try one note first, use a shorter simulation.

```sh
python src/render_trumpet_bach5c.py --only Bb4 --probe 0.45
python src/render_trombone_detailed_scale.py --only Bb3 --probe 0.45
python src/render_tuba_scale.py --only Bb1 --probe 0.45
```

The probe commands save diagnostics for the selected notes. The full-scale commands produce the assembled WAVs.

## Run a CPU model

For the linear bore coupled to the nonlinear lips, generate the saved performance and its diagnostic plot.

```sh
python src/render_scale.py --render
python src/verify.py
```

For the one-dimensional waveguide, including nonlinear return-wave feedback, run the following.

```sh
python src/render_bidirectional.py
python src/verify_bidirectional.py
python src/verify_bidirectional.py --full
python src/plot_bidirectional.py
```

The prescribed-wave approximation is available through `python src/render_nonlinear.py` and `python src/verify_nonlinear.py --full`. Its lip excitation is calculated with the linear resonator and then propagated through the nonlinear model. The bidirectional and axisymmetric renderers calculate the returning wave within their evolving acoustic state.

## Models and geometry

The York trumpet bore and Conn 8H taper are digitized from Campbell, Gilbert, and Myers, *The Science of Brass Instruments*. The trumpet uses a published Bach 5C contour. The trombone's 5AL profile is reconstructed from a cup drawing and calibrated to length, volume, and detached-mouthpiece resonance. The tuba combines an Alexander C-tuba bore scaled to BB♭, Alexander 163 nominal valve and bell dimensions, and a BB-67C4 mouthpiece reconstructed from dimensional and volume measurements. [Geometry and sources](docs/references.md)

The CUDA solver evolves density, axial and radial momentum, and total energy using annular finite volumes, MUSCL reconstruction, an HLLC flux with HLLE fallback, and two-stage SSP Runge–Kutta integration. This conservative formulation can capture shock waves. Lip flow and cup pressure are coupled through the instantaneous returning acoustic characteristic.

The three-octave renderers use local circular-duct thermoviscous losses and four aperture modes at the bell. An equivalent sphere supplies the reflected wave and the radiated pressure at 1 m on axis. The trumpet also includes a cubic transverse lip stiffness. The unfolded bore and reduced exterior preserve axisymmetry. Curved valve passages, wall vibration, resolved turbulent jets, the player's vocal tract, and room acoustics require additional models. [Equations and assumptions](docs/models.md)

## Numerical checks

The shared Euler solver includes a Sod shock-tube check, conservation checks, and a cylindrical radial-mode check. The exterior validator checks stability, passivity, reciprocity, and agreement between two implementations of the spherical-wave reflection.

```sh
python src/verify_axisymmetric.py
python src/verify_bell_exterior.py
```

Checks on the rendered WAVs put the late-sustain pitch error within about 12 cents for the Bach 5C trumpet, within 12 cents for all three 5AL trombone variants, and within 10 cents for the tuba. The reports also identify weak cycle-to-cycle variations where present. These checks describe the rendered performances, whose geometries and player settings have been reconstructed or fitted. [Results and numerical checks](docs/results/)

## Repository layout

```text
src/             Simulation, rendering, verification, and plotting code
data/            Bore profiles, mouthpiece data, performance presets, source records
samples/         Selected WAV files and their regeneration instructions
docs/models.md   Model equations, assumptions, and validation scope
docs/references.md
                 Books, papers, and the role of each reference
docs/results/    Saved measurements and verification reports
docs/figures/    Selected diagnostic plots
docs/ja/         Original model notes in Japanese
output/          Generated audio, traces, and diagnostics, created locally
tmp/             Local caches and working files, created locally
```

## References

The main sources are Campbell, Gilbert, and Myers (2021) for instrument acoustics and geometry, Adachi and Sato (1996) for two-coordinate lip motion, and more recent work on artificial lips and thermoviscous bore modeling. The current implementation combines these ingredients with the assumptions described above. Full citations, links, figure numbers, and how each source is used are collected in [docs/references.md](docs/references.md).
