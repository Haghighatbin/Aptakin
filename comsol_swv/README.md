
# Aptakin — COMSOL SWV Frequency-Response Model (Electrochemical Aptamer Sensor)

**Folder:** `comsol-swv/`</br>
**Model file:** `aptamer_swv.mph` (COMSOL >5.6)</br>
**Purpose:** Predict the **square-wave voltammetry (SWV)** current of a **surface-immobilised aptamer–MB** biosensor as a function of interrogation **frequency** and **analyte concentration**, for rational selection of operating conditions and probe kinetics.

## Executive summary

We model a gold disk electrode functionalised with thiolated DNA aptamers bearing a **methylene blue (MB)** redox label. Target binding drives a conformational change $[unfold (open) ↔ fold (close)]$ that modulates the **MB–electrode distance** and thus the **heterogeneous electron-transfer rate**. Instead of explicitly time-resolving SWV pulses, we collapse SWV to a **frequency-aware exchange current** $i_0(f)$ that captures:</br>
> (*i*) missed slow transitions at high (f); </br>
> (*ii*) low-frequency kinetics-limited behaviour;</br> 
> (*iii*) high-frequency attenuation from the **double layer**.</br>

Outputs: $I(f, C)$, $(\Delta I)$ (bound − unbound), and the **optimal frequency window** $(f_{\text{opt}})$.

---

## 1) Physical system and assumptions

* **Electrode:** Au disk, radius $(r_\mathrm{e} = 1~\text{mm})$.
* **Domain:** 2D axisymmetric electrolyte (PBS, 150 mM), outer radius $(r_\mathrm{dom}=5~\text{mm})$.
* **Surface layer:** aptamer monolayer with MB labels, surface density $(\Gamma_{\mathrm{MB}})$.
* **Analyte:** small-molecule therapeutic, bulk $(C \in [0.1,100]~\mu\text{M})$.
* **Electrolyte physics:** **Secondary Current Distribution** (SCD). High ionic strength ⇒ migration dominates diffusion; no species consumption at the electrode.
* **Temperature:** $(T=298.15~\text{K})$.
* **SWV handling:** no explicit pulse train; use a frequency surrogate embedded in $(i_0(f))$.
* **Reference potential:** formal potential $(E^0)$; interrogation at $(E\approx E^0)$ for differential response.

---

## 2) Mathematical formulation

### 2.1 Two-state surface kinetics (aptamer conformation)

Let $\theta_\mathrm{open}(t)$ and $\theta_\mathrm{closed}(t)=1-\theta_\mathrm{open}(t)$ be the fractional coverages. Surface switching driven by target concentration $(C)$ (bulk, assumed uniform due to $SCD$):

$$
\frac{d\theta_\mathrm{open}}{dt} = k_\mathrm{off},(1-\theta_\mathrm{open}) - k_\mathrm{on},C,\theta_\mathrm{open}
$$

Steady state:
$$
\theta_\mathrm{open}^\ast = \frac{k_\mathrm{off}}{k_\mathrm{off}+k_\mathrm{on}C},\qquad
\theta_\mathrm{closed}^\ast = \frac{k_\mathrm{on}C}{k_\mathrm{off}+k_\mathrm{on}C}
$$

**Parameters used (baseline):** $k_\mathrm{on}=10^4~\text{M}^{-1}\text{s}^{-1}$, $k_\mathrm{off}=0.02~\text{s}^{-1}$

> In the “Route A” stationary workflow we use $(\theta^\ast)$. In “Route B” we time-integrate to steady state and pass values to the stationary step.

---

### 2.2 Distance-dependent electron transfer (Methylene Blue [MB] label)

Electron transfer scales exponentially with donor–electrode distance (d):

$$
k_\mathrm{et}(d) = k_\mathrm{et}^0 ,\exp[-\beta,(d-d_0)]$$

with tunnelling decay $(\beta \approx 1.1~\text{\AA}^{-1})$,
$(d_\mathrm{open} \approx 6.5~\text{nm})$, $(d_\mathrm{closed} \approx 3.2~\text{nm})$, reference $(d_0 = d_\mathrm{closed})$.</br>

State-weighted ET rate:<br/>

$$
\bar{k}*\mathrm{et} = \theta*\mathrm{open},k_\mathrm{et}(d_\mathrm{open}) + \theta_\mathrm{closed},k_\mathrm{et}(d_\mathrm{closed})
$$

---

### 2.3 SWV frequency surrogate (missed transitions at high (f))

Define $(\omega=2\pi f)$. An interrogation that alternates too quickly relative to ET cannot fully relax; we approximate the **effective** rate as:
$$
k_\mathrm{eff}(f) = \frac{\bar{k}*\mathrm{et},\omega}{\omega + \bar{k}*\mathrm{et}}
= \frac{1}{\frac{1}{\bar{k}*\mathrm{et}} + \frac{1}{\omega}}
$$</br>
the harmonic mean of $(\bar{k}*\mathrm{et})$ and $(\omega)$.
Limits: $(k_\mathrm{eff}\to \bar{k}*\mathrm{et})$ for $(\omega\gg\bar{k}*\mathrm{et})$ (fast interrogation does **not** speed ET), and $(k_\mathrm{eff}\to \omega)$ for $(\omega\ll\bar{k}_\mathrm{et})$ (interrogation limits).

---

### 2.4 Double-layer attenuation (RC low-pass)

At high (f), capacitive charging suppresses the Faradaic component. Using $(R_\mathrm{ct} = \dfrac{RT}{nF,i_0^\mathrm{raw}})$ (local linearisation) and surface capacitance $(C_\mathrm{dl})$, we apply a first-order low-pass magnitude:


$$A(\omega) = \frac{1}{\sqrt{1+(\omega R_\mathrm{ct} C_\mathrm{dl})^2}}$$


---

### 2.5 Exchange current and interfacial kinetics (Butler–Volmer)

For a surface-confined redox with surface density $(\Gamma_{\mathrm{MB}})$ and electron number $(n=1)$:

$$
i_0^\mathrm{raw}(f) = n F \Gamma_{\mathrm{MB}},k_\mathrm{eff}(f),\qquad
i_0(f) = i_0^\mathrm{raw}(f),A(\omega).
$$

The interfacial current density:
$$
i_\mathrm{BV} = i_0(f)\left[\exp!\left(\frac{\alpha_a nF\eta}{RT}\right) - \exp!\left(-\frac{\alpha_c nF\eta}{RT}\right)\right],\quad
\eta = ( \phi_s - \phi_l ) - E^0.
$$

We evaluate near $(E\approx E^0) (small (|\eta|))$ to emulate differential SWV response.

---

### 2.6 Total current (2D axisymmetric)

The electrode current is the **axisymmetric line integral** over the disk edge:

$$
I = \int_{\partial \Omega_\mathrm{elec}} i_n . 2\pi r . \mathrm{d}r
$$

where $(i_n)$ is the normal current density (A/m²) from SCD.

---

## 3) Parameters (baseline values)

| Symbol                 |             Value | Units       | Note                       |
| ---------------------- | ----------------: | ----------- | -------------------------- |
| $(F)$                    |             96485 | C·mol⁻¹     | Faraday                    |
| $(R)$                    |             8.314 | J·mol⁻¹·K⁻¹ | Gas constant               |
| $(T)$                    |            298.15 | K           |                            |
| $(n)$                    |                 1 | –           | MB redox                   |
| $(\Gamma_{\mathrm{MB}})$ | $(1\times10^{-11})$ | mol·cm⁻²    | → $(1\times10^{-7})$ mol·m⁻² |
| $(k_\mathrm{on})$        |            (10^4) | M⁻¹·s⁻¹     |                            |
| $(k_\mathrm{off})$       |              0.02 | s⁻¹         |                            |
| $(\beta)$                |             (1.1) | Å⁻¹         | $(=1.1\times10^{10})$ m⁻¹    |
| $(d_\mathrm{open})$      |               6.5 | nm          |                            |
| $(d_\mathrm{closed})$    |               3.2 | nm          |                            |
| $(k_\mathrm{et}^0)$      |               500 | s⁻¹         | pre-exponential (tunable)  |
| $(C_\mathrm{dl})$        |                20 | μF·cm⁻²     | → (0.02) F·m⁻²             |
| $(\alpha_a,\alpha_c)$    |          0.5, 0.5 | –           |                            |
| $(E^0)$                  |             −0.26 | V (Ag/AgCl) | set per dataset            |
| $(\sigma_\mathrm{elec})$ |               1.5 | S·m⁻¹       | PBS (150 mM)               |
| $(r_\mathrm{e})$         |                 1 | mm          | electrode radius           |
| $(r_\mathrm{dom})$       |                 5 | mm          | domain radius              |

**Sweep ranges:** $(f \in [10,1000]~\text{Hz})$ (log), $(C \in [0.1,100]~\mu\text{M})$ (log).

---

## 4) COMSOL implementation (concise HOWTO)

### 4.1 Geometry (2D Axisymmetric)

* Rectangle: width $(r_\mathrm{dom})$, height $(r_\mathrm{dom})$.
* Bottom edge $(z=0)$: segment $(0\le r\le r_\mathrm{e})$ is the **electrode**.

### 4.2 Physics

* **Secondary Current Distribution (SCD):** electrolyte domain with conductivity $(\sigma_\mathrm{elec})$.
* **Electrode boundary:** **Butler–Volmer** with **Exchange current density** $(i_0(f))$, **$(\alpha_a,\alpha_c)$**, and **Double-layer capacitance** $(C_\mathrm{dl})$.
* **Applied potential:** set solid potential $(\phi_s = E)$ (use $(E\approx E^0)$ for frequency scans).
* All other boundaries: insulation/symmetry; outer boundary: fixed bulk (C) is handled parametrically (no transport PDEs in SCD).

### 4.3 Variables (Global → Variables)

Define:

* $(\theta^\ast_\mathrm{open}, \theta^\ast_\mathrm{closed})$ using $(k_\mathrm{on},k_\mathrm{off},C)$
* $(k_\mathrm{et}(d_\mathrm{open/closed}))$, $(\bar{k}_\mathrm{et})$
* $(\omega=2\pi f)$, $(k_\mathrm{eff}(f))$, $(i_0^\mathrm{raw})$, $(R_\mathrm{ct})$, $(A(\omega))$, $(i_0(f))$

*(If using the transient pre-step, replace $(\theta^\ast)$ by the time-dependent $(\theta)$ from **Global ODEs and DAEs**.)*

### 4.4 Mesh

Physics-controlled is sufficient. Optional boundary layer on the electrode: 4–6 layers, first thickness $(5\times10^{-10})$ m.

### 4.5 Studies

* **Route A (recommended for sweeps):** **Stationary** (SCD only) with a **Parametric Sweep** over (f) and (C).
* **Route B:** **Time Dependent** (global ODE for $(\theta)$) to steady state → **Stationary** SCD using those values.

### 4.6 Derived values and plots

* **Line Integration** over electrode boundary with **Use axisymmetry** to compute total current (I).
* **1D Plot Group → Global**: (I) vs (f) (and separate curves for (C)).
* Export tables/figures to `exports/`.

---

## 5) Automation with LiveLink (MATLAB)

**File:** `run_swv_grid.m`

* Loads `aptamer_swv.mph`.
* Sets parameters ((f, C)) on a log grid (25 × 15 = 375 sims).
* Solves the Stationary study (`std1`).
* Evaluates total current with your integration operator, e.g. `-intI(j_s)` or a **Line Integration** expression.
* Saves `results/swv_grid.mat` and `exports/heatmap.png`.

**Run:**

```matlab
cd comsol-swv
run_swv_grid
```

> We'll have to Ensure the expression in `currentExpr` matches what works in **Results → Global Evaluation** inside COMSOL.

---

## 6) Validation protocol (recommended)

1. **Peak current vs frequency:** At a fixed $(C)$ (e.g., 10 µM), compare simulated $(I(f))$ against experimental SWV peak currents; report $(R^2)$.
2. **Dose–response at $(f_{\text{opt}})$:** Extract $(I(C))$, fit with $(\theta^\ast_\mathrm{closed}(C))$ and report KD error (%).
3. **Kinetic extraction:** From a small **Time Dependent** step in (C), fit relaxation rate $(k_\mathrm{off}+k_\mathrm{on}C)$; compare with orthogonal measurements (e.g., SPR).

**Deliverables:** `validation/compare_sim_vs_exp.ipynb` prints $(R^2)$ and KD error; figures saved to `docs/figs/`.

---

## 7) Sensitivity analysis

Vary each parameter ±50% (one-at-a-time or Latin Hypercube) and map the effect on:

* $(f_{\text{opt}})$ shift $(\Delta \log_{10} f)$,
* $\max \Delta I$ (contrast),
* KD bias from dose–response fitting.

Key levers typically: $k_\mathrm{et}^0$, $\beta$, $C_\mathrm{dl}$, $\Gamma_{\mathrm{MB}}$, $k_\mathrm{on}/k_\mathrm{off}$.

---

## 8) Limitations and extensions

* **SWV surrogate:** captures major trends but omits pulse-by-pulse dynamics and harmonic content. For spot checks, implement explicit SWV with **Events** (step-hold waveform).
* **Mass transport:** SCD neglects concentration gradients of the analyte; appropriate for small-molecule targets at high supporting electrolyte and low consumption. If needed, add **Transport of Diluted Species** with surface reaction sink/source.
* **Potential dependence of conformation:** If electrostatic gating is relevant, extend $(\theta)$ rates to include (E) dependence.
* **Noise & instrumentation:** Add a readout noise model or EIS-derived $C_\mathrm{dl}(E)$ for instrument-specific predictions.

---

## 9) Reproduction checklist

* COMSOL 6.1 + Electrochemistry Module; (optional) LiveLink for MATLAB.
* Open `aptamer_swv.mph` → run **Study 1** → plot **I vs f**.
* MATLAB: `run_swv_grid` → inspect `results/swv_grid.mat` and `exports/heatmap.png`.
* (Optional) Add experimental CSVs to `validation/experimental_data` → run notebook.

---

## 10) File map

```
comsol-swv/
├─ aptamer_swv.mph          # COMSOL model (Git LFS)
├─ run_swv_grid.m           # LiveLink sweep script
├─ params/                  # (optional) YAML/JSON/TXT parameter sets -> params, aptamer_var
├─ results/                 # MAT/CSV (Git LFS)
├─ images/                  # Sample images (PNG)
├─ exports/                 # PNG/PDF reports (Git LFS)
└─ README.md                # this document
```

---

## 11) Symbols and units (quick reference)

| Symbol | Meaning | Units | Baseline / Notes |
|--------|---------|-------|------------------|
| $F$ | Faraday constant | $C·mol^{-1}$ | $96485$ |
| $R$ | Gas constant | $J·mol^{-1}·K^{-1}$ | $8.314$ |
| $T$ | Temperature | $K$ | $298.15$ |
| $n$ | Electrons transferred | – | $1 (MB-label)$ |
| $\Gamma_{\mathrm{MB}}$ | Surface density of MB redox sites | $mol·m^{-2}$ | $1\times10^{-7} (from 1\times10^{-11} mol·cm^{-2}$) |
| $k_{\mathrm{on}}$ | Association rate | $M^{-1}·s^{-1}$ | $1\times10^{4}$ |
| $k_{\mathrm{off}}$ | Dissociation rate | $s^{-1}$ | $0.02$ |
| $\theta_{\mathrm{open}}, \theta_{\mathrm{closed}}$ | Surface state fractions | – | $\theta_{\mathrm{closed}}=1-\theta_{\mathrm{open}}$ |
| $\beta$ | Tunnelling decay constant | $m^{-1}$ | $1.1\times10^{10}$ ($1.1$ $Å^{-1}$) |
| $d_{\mathrm{open}}, d_{\mathrm{closed}}$ | Donor–electrode distances | m | $6.5 nm, 3.2 nm$ |
| $d_{0}$ | Reference distance | m | $d_{\mathrm{closed}}$ |
| $k_{\mathrm{et}}^{0}$ | ET pre-exponential | $s^{-1}$ | $500 (tunable)$ |
| $k_{\mathrm{et}}(d)$ | Distance-dependent ET rate | $s^{-1}$ | $k_{\mathrm{et}}^{0}\exp[-\beta(d-d_{0})]$ |
| $\bar{k}_{\mathrm{et}}$ | State-weighted ET rate | $s^{-1}$ | $\theta_{\mathrm{open}}k_{\mathrm{et}}(d_{\mathrm{open}})+\theta_{\mathrm{closed}}k_{\mathrm{et}}(d_{\mathrm{closed}})$ |
| $f$ | Interrogation frequency | $Hz$ | $sweep: 10–1000$ |
| $\omega$ | Angular frequency | $s^{-1}$ | $2\pi f$ |
| $k_{\mathrm{eff}}(f)$ | Effective ET under SWV | $s^{-1}$ | $\dfrac{\bar{k}_{\mathrm{et}}\omega}{\omega+\bar{k}_{\mathrm{et}}}$ |
| $C_{\mathrm{dl}}$ | Double-layer capacitance (surface) | $F·m^{-2}$ | $0.02 (20 μF·cm^{-2}$) |
| $R_{\mathrm{ct}}$ | Charge-transfer resistance (surface) | $Ω·m^{2}$ | $\dfrac{RT}{nF\,i_{0}^{\mathrm{raw}}}$ |
| $A(\omega)$ | RC attenuation magnitude | – | $\dfrac{1}{\sqrt{1+(\omega R_{\mathrm{ct}} C_{\mathrm{dl}})^{2}}}$ |
| $i_{0}^{\mathrm{raw}}(f)$ | Raw exchange current density | $A·m^{-2}$ | $nF\Gamma_{\mathrm{MB}}\,k_{\mathrm{eff}}(f)$ |
| $i_{0}(f)$ | Attenuated exchange current | $A·m^{-2}$ | $i_{0}^{\mathrm{raw}}(f)\,A(\omega)$ |
| $\alpha_{a}, \alpha_{c}$ | Transfer coefficients | – | $0.5, 0.5$ |
| $\eta$ | Overpotential | $V$ | $\eta=(\phi_{s}-\phi_{l})-E^{0}$ |
| $E^{0}$ | Formal potential | $V$ | $−0.26$ *vs.* $Ag/AgCl$ (example) |
| $i_{\mathrm{BV}}$ | Butler–Volmer current density | $A·m^{-2}$ | see equation below |
| $\sigma_{\mathrm{elec}}$ | Electrolyte conductivity | $S·m^{-1}$ | $1.5$ ($PBS$ $150 mM$) |
| $r_{\mathrm{e}}$ | Electrode radius | $m$ | $1 mm$ |
| $r_{\mathrm{dom}}$ | Domain radius | $m$ | $5 mm$ |
| $I$ | Total electrode current | $A$ | Axisymmetric line integral |

### Core Equations

**Aptamer kinetics (surface, two-state):**

$$
\frac{d\theta_{\mathrm{open}}}{dt}=k_{\mathrm{off}}(1-\theta_{\mathrm{open}})-k_{\mathrm{on}}C\,\theta_{\mathrm{open}},\quad
\theta_{\mathrm{closed}}^{\ast}=\frac{k_{\mathrm{on}}C}{k_{\mathrm{on}}C+k_{\mathrm{off}}}
$$

**Distance-dependent ET & SWV surrogate:**

$$
k_{\mathrm{et}}(d)=k_{\mathrm{et}}^{0}e^{-\beta(d-d_{0})},\quad
\bar{k}_{\mathrm{et}}=\theta_{\mathrm{open}}k_{\mathrm{et}}(d_{\mathrm{open}})+\theta_{\mathrm{closed}}k_{\mathrm{et}}(d_{\mathrm{closed}}),\quad
k_{\mathrm{eff}}(f)=\frac{\bar{k}_{\mathrm{et}}\omega}{\omega+\bar{k}_{\mathrm{et}}}
$$

**Double-layer attenuation & exchange current:**

$$
A(\omega)=\frac{1}{\sqrt{1+(\omega R_{\mathrm{ct}}C_{\mathrm{dl}})^{2}}},\qquad
i_{0}(f)=nF\Gamma_{\mathrm{MB}}k_{\mathrm{eff}}(f)\,A(\omega)
$$

**Butler–Volmer kinetics (evaluated near $E\approx E^{0}$):**

$$
i_{\mathrm{BV}}=i_{0}(f)\left[\exp\!\left(\frac{\alpha_{a}nF\eta}{RT}\right)-\exp\!\left(-\frac{\alpha_{c}nF\eta}{RT}\right)\right],\quad
\eta=(\phi_{s}-\phi_{l})-E^{0}
$$

**Total current (2D axisymmetric):**

$$
I=\int_{\partial\Omega_{\mathrm{elec}}} i_{n}\,2\pi r\,\mathrm{d}r
$$

where $i_{n}$ is the normal current density from SCD.

---

### Acknowledgements

This COMSOL module complements the Python side of **Aptakin** and is intended for transparent, reproducible, and interview-friendly demonstration of continuum electrochemistry modelling with explicit scientific rationale.

---
### Developer
**Amin Haghighatbin**<br/>
email: aminhb@tutanota.com
