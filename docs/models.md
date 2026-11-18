# Physical models

The renderers use SI units internally. Sounding pitches are used throughout. Each scale is assembled from separate notes with a prescribed blowing-pressure envelope and a fixed bore geometry during the note.

## Lip motion and the entrance boundary

`src/trumpet.py` and `src/axisymmetric_euler.cu` implement Adachi and Sato's two-coordinate lip model. The upper lip tip has axial and transverse coordinates `x` and `y`. Pressure forces drive the motion, while a spring and damping oppose displacement and velocity. Contact adds a restoring force and stronger damping when `y < 0`.

The open gap area and swept-volume flow are

```text
S_gap = max(2*b*y, 0)
U_swept = b * (x*dy/dt - (y-y_joint)*dx/dt)
U_total = U_air + U_swept
```

The gap flow includes fluid inertia and a signed quadratic pressure loss. Given the returning wave and lip state, the implementation solves the resulting quadratic flow equation at the entrance.

For the axisymmetric solver, the area-weighted returning characteristic and entrance pressure satisfy

```text
h = area_mean(p_prime - rho0*c0*u_z)
p_cup = h + (rho0*c0/S_cup)*U_total
```

The tube supplies `h` at each fluid stage. The cup and gap pressures act on the subsequent lip update. The entrance uses a small-Mach-number acoustic characteristic approximation even though the bore interior evolves with the nonlinear Euler equations.

The default trumpet lip width is 7 mm, lip thickness 2 mm, pivot height 4 mm, equilibrium axial coordinate 1 mm, and cup entrance area 230 mm². The mass and spring depend on the lip-frequency control using Adachi and Sato's parameterization. The trombone uses larger effective dimensions and a different mass parameter, listed in `trombone_axisymmetric.LIPS`. Per-note blowing pressures, equilibrium openings, and lip frequencies are fitted performance settings.

## Cubic stiffness in the current trumpet

The transverse restoring acceleration in `src/axisymmetric_euler.cu` is

```text
eta = y - y_equilibrium
a_spring_y = -omega_y^2 * eta * (1 + alpha*eta^2)
```

`src/render_trumpet_modern.py` sets `alpha = 8e4 m^-2`. At a displacement of 1 mm, the cubic restoring term is 8% of the linear term.

This extension is informed by Fréour and colleagues' reduced artificial-lip model with cubic stiffness and by Inomata and colleagues' study of geometric stiffness under blowing pressure. Combining this term with the Adachi–Sato two-coordinate system is a modeling choice. Its coefficient has not been identified from measurements of a human player's lips. The original contact and flow models are retained.

## Instrument geometry

The trumpet body is reconstructed from the York B-flat trumpet, Grand Rapids, 1925–27, EU(5772), in Figure 7.52 of *The Science of Brass Instruments*. It is sampled from a raster graph. The approximate reading resolution is 2.8 mm axially and 0.43 mm in diameter, with greater uncertainty at crossings and near the bell rim.

The trumpet mouthpiece is reconstructed as a rounded converging cup, a 3.7 mm diameter throat 6 mm long, and a 55 mm backbore. Its cup depth was calibrated to the book's 806 Hz detached-mouthpiece resonance using the original linear acoustic model. This single resonance constrains the reconstruction only partially. The newer axisymmetric and loss models need not predict exactly that calibration frequency.

Valve tubing is represented by lengths of approximately cylindrical bore corresponding to the selected semitone extensions. Per-note tuning offsets adjust the effective cylindrical length. A negative total offset removes a portion of the tuning section.

The trombone uses the Conn 8H profile in Figure 4.22, the cup profile in Figure 4.32b, and representative Denis Wick 5AL dimensions. The reconstruction constrains mouthpiece length to 80 mm, volume to 12 mL, and a detached resonance to 535 Hz in the linear calibration. Slide movement extends both legs of the unfolded bore. Per-note `trim_m` provides an additional effective-length correction.

The instrument centerline is unfolded into an axisymmetric profile. The mesh therefore represents a radius that varies with distance along the bore. The valve ports and U-bends would require additional geometry to resolve their three-dimensional effects.

## Nonlinear axisymmetric flow

The CUDA solver in `src/axisymmetric_euler.cu` evolves density, axial momentum, radial momentum, and total energy of an ideal gas with `gamma = 1.4`, `rho0 = 1.2 kg/m³`, and `c0 = 343 m/s`.

```text
Q = (rho, rho*u_z, rho*u_r, E)
F_z = (rho*u_z, rho*u_z^2 + p, rho*u_z*u_r, (E+p)*u_z)
F_r = (rho*u_r, rho*u_z*u_r, rho*u_r^2 + p, (E+p)*u_r)
p = (gamma-1) * (E - rho*(u_z^2+u_r^2)/2)

d(r*Q)/dt + d(r*F_z)/dz + d(r*F_r)/dr = (0, 0, p, 0)
```

The mesh uses annular finite volumes with conical side faces. The solver uses MUSCL reconstruction with an MC limiter, an HLLC Riemann flux, an HLLE fallback for inadmissible intermediate states, and SSP-RK2 integration. The geometric pressure source balances the pressure flux at rest. Wave-speed and positivity checks stop a run when its state exceeds the configured limits.

This conservative formulation allows steep compression waves and shock waves. A Sod shock-tube test checks the flux and conservation. The mesh does not resolve the physical microscopic thickness of an acoustic shock.

The trumpet's default axial spacing is 3.5 mm, with finer spacing in the mouthpiece and final flare. The bore has 2–16 radial rings in the nominal open-fingering mesh. The default trombone spacings are 6 mm axially and a nominal 4.5 mm radially, also with local refinement. The global time step is determined by a CFL estimate with a wave-speed reserve. Boundary histories are recorded at 192 kHz and radiated audio is filtered and resampled to 48 kHz.

## Thermoviscous closure

The earlier axisymmetric presets use five causal relaxation terms to approximate a loss proportional to the square root of frequency.

The current trumpet uses local circular-duct viscous and thermal responses following the Bessel-function treatment also used by Jackson and Scavone. With Laplace variable `s`, local radius `a`, kinematic viscosity `nu`, and thermal diffusivity `D`,

```text
w_v = s*a^2/nu
w_T = s*a^2/D
F(w) = 2*I_1(sqrt(w)) / (sqrt(w)*I_0(sqrt(w)))
D_v(s) = s*F(w_v)/(1-F(w_v))
D_T(s) = (gamma-1)*s*F(w_T)
```

`I_n` denotes a modified Bessel function. `src/thermoviscous.py` fits the dynamic parts with ten positive relaxation weights and positive pole rates for each response. The viscous DC term `8*nu/a^2` is included separately. The parameters are `nu = 1.54e-5 m²/s` and `D = 2.17e-5 m²/s`. Maximum relative errors on the fitting band are approximately 1.56% for the viscous dynamic response and 2.31% for the thermal response.

The resulting momentum and pressure relaxation is distributed across each section. It is a local linear circular-duct closure applied to a nonlinear fluid calculation. Actual thermal and viscous boundary layers are represented through this closure rather than spatially resolved at the wall.

## Bell reflection, inflow, and audible output

The original axisymmetric presets use a passive first-order RL load with a low-frequency end correction.

The current trumpet uses Silva and colleagues' causal reflection approximation for an unflanged circular opening. For `x = s*a/c0`,

```text
R(s) = -(1 + 0.167*x) / (1 + 1.393*x + 0.457*x^2)
```

The filter acts on the outgoing mean pressure wave at the bell. Other radial components use a local outgoing characteristic. The current trumpet sets the entropy of gas entering from an external reservoir to the reference value, while exiting gas carries the local interior entropy.

The radiated pressure is estimated at a distance `d = 1 m` as

```text
p_radiated(t) = rho0/(4*pi*d) * dU_bell(t)/dt
```

This compact-monopole estimate omits an external radiation grid and angular directivity. Soares and colleagues' 2025 work provides an example of an axisymmetric modal method that includes an exterior field and directivity. Their reduced-order exterior model is a reference for further development.

## Earlier CPU models

| Renderer | Acoustic calculation | Lip feedback |
| --- | --- | --- |
| `render_scale.py` | Lossy, linear transfer matrices converted to a causal reflection response | Linear returning acoustic response |
| `render_nonlinear.py` | Generalized Burgers propagation of a prescribed incident wave | Excitation supplied by the linear resonator |
| `render_bidirectional.py` | Area-junction waveguide, causal wall losses, and Burgers self-distortion in each direction | Actual evolving returning wave |

The bidirectional CPU model neglects quadratic mixing between the two counterpropagating waves. Its assumptions and source records are in `data/sources_bidirectional.json`. The prescribed-wave model's scope is recorded in `data/nonlinear_sources.json`.

## Validation and sample results

The committed JSON files in `docs/results/` record the published samples and selected numerical checks. Generated results from a new run are written under `output/`.

- `verify_axisymmetric.py` checks a resting state, mass and energy conservation, the Sod shock tube at two resolutions, and a cylindrical radial eigenmode.
- `verify_trumpet_modern.py` checks the bell reflection magnitude, positive thermoviscous fits, weak-pulse input impedance, and a B♭4 run with return-wave coupling disabled.
- `verify_trombone.py` checks mouthpiece geometry and resonance, slide lengths, and bore impedance. Its `--coupling` option compares B♭3 with disabled feedback and a refined mesh.
- The CPU verification scripts provide analytical, refinement, and audio measurements appropriate to their respective acoustic approximations.

In the current trumpet sample, the maximum late-sustain pitch error is 8.2 cents. D4, E♭4, and F4 show cycle-to-cycle variation. In the tested B♭4 configuration, the cup pressure is about 5.37 kPa RMS and sustained oscillation disappears with return-wave feedback disabled.

Updating the wall and bell models on the same nominal open-fingering mesh shifts the first six weak-signal resonances by less than 0.7 Hz. The first impedance peak decreases from about 113.8 to 96.8 MPa·s/m³. These are comparisons between numerical models. Mesh sensitivity, graph-derived dimensions, and the reconstructed mouthpieces must be considered before drawing conclusions about a specific physical instrument.

The model and its coefficients can be explored independently of the saved performance presets. Matching a selected musical scale does not determine all physical parameters or establish agreement with every playing regime. The cubic stiffness, local boundary-layer closure, and compact radiation approximation require further validation against measured lips, impedances, and radiated sound.
