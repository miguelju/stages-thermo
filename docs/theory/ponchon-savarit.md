# Ponchon–Savarit — the enthalpy–composition method

![Ponchon–Savarit construction on the enthalpy–composition (H–x–y) diagram for a
benzene–toluene 95/5 split at 101.325 kPa: blue saturated-liquid and red
saturated-vapor enthalpy curves, green tie lines (equilibrium stages), purple
pole (operating) lines, and the two difference points Δ_D (top) and Δ_B (bottom)
joined through the feed point F by the dotted balance line.](img/ponchon-savarit.png)

Rung 2 of the pedagogical ladder. Where [McCabe–Thiele](mccabe-thiele.md) assumes
**constant molal overflow** (CMO) and steps a staircase on an $x$–$y$ square,
Ponchon–Savarit closes the **energy balance exactly** on the
enthalpy–composition (H–x–y) diagram. It is the method you need when the molar
latent heats are unequal or the heat of mixing is large — most dramatically for
ammonia–water.

**Method and references.** Ponchon, M. *La Technique Moderne* **1921**, 13, 20
and 55. Savarit, R. *Arts et Métiers* **1922**. Implementation-level treatment:
Seader, Henley & Roper (S&H), *Separation Process Principles*, Ch. 7. See
[`docs/references.md`](../references.md) for the full list.

Code: [`engine/src/binary/ponchon_savarit.rs`](../../engine/src/binary/ponchon_savarit.rs),
built on the enthalpy-augmented curve in
[`engine/src/binary/equilibrium.rs`](../../engine/src/binary/equilibrium.rs)
(`EnthalpyCurve`) and the per-phase molar enthalpies in
[`engine/src/thermo.rs`](../../engine/src/thermo.rs) (`ThermoSystem::phase_enthalpy`).

## Units

Enthalpies in **kJ/kmol**, compositions are light-component mole fractions
(dimensionless), pressure **kPa** (absolute), temperature **K**. Duties are
reported **per mole of feed** (kJ/kmol of feed).

## The diagram

Instead of $y^\ast(x)$ on a unit square, the construction lives on a plot of
**molar enthalpy versus composition** with two saturation curves:

- the **saturated-liquid** curve $h_L(x)$ — the molar enthalpy of a boiling
  liquid of composition $x$ at its bubble temperature $T(x)$;
- the **saturated-vapor** curve $H_V(y)$ — the molar enthalpy of a saturated
  vapor of composition $y$.

A liquid and its incipient vapor share a temperature at the bubble point, so the
pair $(x, y^\ast(x))$ is joined by a **tie line** on the diagram: from
$(x, h_L(x))$ to $(y^\ast(x), H_V(y^\ast(x)))$, both evaluated at $T(x)$. In the
code these arrays are built once by `EnthalpyCurve::from_thermo` (compute them
from a model) or `EnthalpyCurve::from_points` (feed reference data).

## Difference points (poles)

The energy balance is carried by two **difference points** — the poles $\Delta_D$
(top) and $\Delta_B$ (bottom). Writing $Q_C$ for the condenser duty and $Q_R$
for the reboiler duty, with distillate $D$ and bottoms $B$:

$$\Delta_D = \left(x_D,\ h_D + \frac{Q_C}{D}\right), \qquad \Delta_B = \left(x_B,\ h_B - \frac{Q_R}{B}\right).$$

Write the ordinate of $\Delta_D$ as $Q_D^{\prime} = h_D + Q_C/D$. The **reflux
ratio sets the height of the top pole**: reading $R = L/D$ off the H–x diagram as
a ratio of vertical segments through the condenser (S&H Ch. 7),

$$R = \frac{Q_D^{\prime} - H_{V1}}{H_{V1} - h_{L0}} \qquad\Longrightarrow\qquad Q_D^{\prime} = H_{V1} + R (H_{V1} - h_{L0}),$$

where $H_{V1}$ is the saturated-vapor enthalpy at $y = x_D$ (the vapor entering a
total condenser) and $h_{L0}$ is the saturated-liquid reflux enthalpy at $x_D$.

**The overall balance makes the three points collinear.** A material and energy
balance around the whole column places the feed point $F = (z_F, h_F)$ on the
straight line joining $\Delta_D$ and $\Delta_B$. Given $\Delta_D$ and $F$, this
fixes $\Delta_B$: extrapolate the line $\Delta_D$–$F$ to $x = x_B$.

The feed enthalpy follows from the thermal condition $q$ (1 = saturated liquid,
0 = saturated vapor), read off the same two curves:

$$h_F = q \cdot h_L(z_F) + (1 - q) \cdot H_V(z_F).$$

## Duties

With a total condenser the distillate is saturated liquid, $h_D = h_L(x_D)$, and
the lever rule splits one mole of feed into $D/F = (z_F - x_B)/(x_D - x_B)$ and
$B/F = (x_D - z_F)/(x_D - x_B)$. The duties per mole of feed are the pole
ordinates measured from the product enthalpies:

$$\frac{Q_C}{F} = \frac{D}{F} (Q_D^{\prime} - h_D), \qquad \frac{Q_R}{F} = \frac{B}{F} (h_B - Q_B^{\prime}),$$

and the whole-column energy balance closes as

$$F \cdot h_F + Q_R = D \cdot h_D + B \cdot h_B + Q_C.$$

This closure is asserted as an invariant on every converged construction (see the
tests in `ponchon_savarit.rs` and `engine/tests/m2_ponchon_savarit.rs`).

## Stepping (code ↔ diagram)

The staircase mirrors McCabe–Thiele's exactly, with **two substitutions**:

1. **Equilibrium step → tie line.** The liquid $x_n$ leaving stage $n$ is in
   equilibrium with the vapor $y_n$ leaving it: $x_n = x^\ast(y_n)$. This is the
   *same* inverse the McCabe–Thiele stepper uses (`EquilibriumCurve::x_of_y`); on
   the H–x–y diagram it is the tie line through $(x_n, h_L(x_n))$ and
   $(y_n, H_V(y_n))$.

2. **Operating step → pole line.** The passing streams $L_n$ (liquid, $x_n$) and
   $V_{n+1}$ (vapor from below) plus the section's pole $\Delta$ are **collinear**.
   So $y_{n+1}$ is the vapor composition where the straight line through $\Delta$
   and $(x_n, h_L(x_n))$ cuts the saturated-vapor curve. Collinearity of
   $\Delta = (x_p, h_p)$, the liquid point $L = (x_l, h_l)$ and a candidate vapor
   point $V = (y, H_V(y))$ is the 2-D cross product

   $$g(y) = (x_l - x_p)(H_V(y) - h_p) - (h_l - h_p)(y - x_p) = 0.$$

   $H_V$ is piecewise-linear on the sampled vapor grid, so $g$ is piecewise-linear
   and the root inside a bracketing segment is exact — the same sign-change scan
   McCabe–Thiele uses for the q-line/curve crossing.

Stages are stepped **top-down** (stage 1 at the top). The rectifying section uses
the pole $\Delta_D$; the stripping section uses $\Delta_B$. The **feed stage** is
the first stage stepped on $\Delta_B$; the switch happens when the stage liquid
passes the feed composition $z_F$ — which, for a saturated-liquid feed,
coincides with McCabe–Thiele's optimal-feed rule. A partial reboiler is an
equilibrium stage and is counted; v1 supports a **total condenser** (a partial
condenser adds an equilibrium stage above the top tray that is not yet modelled).

## Symbol map (code ↔ textbook)

```text
spec.x_distillate   x_D              distillate light-component fraction
spec.x_bottoms      x_B              bottoms light-component fraction
spec.z_feed         z_F              feed light-component fraction
spec.q              q                feed thermal condition
spec.reflux         R = L/D          external reflux ratio
delta_d             Δ_D = (x_D, Q'_D)  top difference point
delta_b             Δ_B = (x_B, Q'_B)  bottom difference point
q_condenser         Q_C/F            condenser duty per mole feed
q_reboiler          Q_R/F            reboiler duty per mole feed
```

## Consistency with McCabe–Thiele

For a **near-ideal, equal-latent-heat** system (benzene–toluene) CMO is a good
assumption, so the pole lines reproduce the straight CMO operating lines and the
two methods return the same stage count — the M2 consistency check (they agree to
within about a stage; the small residual gap *is* benzene–toluene's real CMO
error, since its latent heats differ by a few percent). Where the heat of mixing
is large the two diverge, and Ponchon–Savarit is the one that is right.

---

## Ammonia–water: how the textbook did it, and the limits of a computed chart

*Motivating the two ammonia–water constructions in
[`notebooks/02-ponchon-savarit.ipynb`](../../notebooks/02-ponchon-savarit.ipynb) —
one computed from a model, one from reference data.*

**The textbook case was done on charts, not on-the-fly thermo.** Ponchon (1921)
and Savarit (1922) devised the enthalpy–composition method as a *graphical*
procedure run on a pre-drawn H–x–y diagram of the binary at hand: nothing is
computed during the construction — you locate the difference points (poles) and
step off stages with a straightedge on a printed chart. For ammonia–water that
chart is the Merkel–Bošnjaković enthalpy–concentration diagram (Merkel &
Bošnjaković, 1929; Bošnjaković, *Technische Thermodynamik*, 1935), built from
**experimental data** — vapor–liquid equilibrium plus **calorimetry** (heats of
solution/mixing, heat capacities, latent heats) — reduced onto one diagram. Its
modern refinements are the correlations still used today: Scatchard et al.
(1947), Macriss et al. (1964), the ASHRAE/EES formulation of Ibrahim & Klein
(1993), and the reference Helmholtz-energy EOS of Tillner-Roth & Friend (1998).
So the "textbook thermo" for ammonia–water *is* a data artifact — a
chart/correlation — and Ponchon–Savarit is the graphical construction performed
on it.

**Why ammonia–water was the showcase.** Its heat of mixing is enormous. Constant
molal overflow — the assumption behind McCabe–Thiele — requires equal molar
latent heats and negligible enthalpy of mixing so that liquid and vapor flows
stay constant along each section. Ammonia–water violates that badly, so
McCabe–Thiele returns the wrong stage count *and* the wrong internal flows.
Ponchon–Savarit closes the energy balance exactly through the chart, so it is the
method you genuinely need here — which is precisely why textbooks introduce it on
ammonia–water.

**Why a computed NRTL γ-φ chart frays.** Regenerating the H–x–y diagram from a
model splits it into pieces, each with its own error: activity coefficients
(NRTL) for the liquid, a cubic EOS for the vapor. The vapor enthalpy is
ideal-gas terms plus an EOS departure; the liquid enthalpy is ideal-gas terms
**minus a per-component condensation (latent-heat) term plus the NRTL excess
enthalpy** — the γ-φ route vle-thermo actually computes. The weak links, in
order: (1) **the vapor at elevated pressure** — the classic charts live at
several bar, where the cubic-EOS vapor enthalpy departure is only approximate;
this dominates the error and is not the activity model's fault; (2) **the
pure-liquid-ammonia reference for a dissolving light gas** — the γ-φ liquid route
reaches the liquid through each pure component's latent heat, a
Clausius–Clapeyron term taken from the slope of the saturation-pressure
correlation, so a few-percent error in that slope lands directly in the liquid
enthalpy; and pure liquid ammonia itself is an awkward datum as temperature
climbs toward its critical point (405.4 K — generator conditions get close);
(3) **the curvature of the excess enthalpy** — the real ammonia–water heat of
mixing is large, asymmetric, and strongly temperature-dependent across the full
range, and a few-parameter NRTL cannot match it at the water-rich and
ammonia-rich ends at once. Net: NRTL gets a few percent at moderate pressure — a
real improvement over Wilson or van Laar — but not a pixel-match to the
Bošnjaković chart.

**The "UNIQUAC is in the good models" trap.** It is tempting to reach for UNIQUAC
because the celebrated ammonia–water model uses it — but that model (Thomsen &
Rasmussen, 1999) is *extended* UNIQUAC: plain UNIQUAC plus a Debye–Hückel term
plus speciation. All the accuracy comes from those additions, not from the
local-composition kernel. On a single binary fit, plain UNIQUAC has only two
adjustable energy parameters (its $r$ and $q$ are fixed structural constants)
against NRTL's three, and its size-asymmetry advantage is wasted on two small
molecules — so plain UNIQUAC is, if anything, slightly *less* flexible than NRTL
here and no closer to the chart.

**What actually would get closer,** in increasing order of effort: (1) feed the
reference correlation's data (Ibrahim–Klein / Tillner-Roth) directly and run the
pure graphical construction on it — matches the textbook chart, needs no new
model; (2) implement *extended* UNIQUAC or a Helmholtz-energy EOS — the models
that genuinely reproduce the chart, but each a large undertaking.

**Why this project does not build the specialized models.** Their distinguishing
capability serves nothing else on the roadmap. Extended UNIQUAC's
electrolyte/speciation machinery is used by exactly one planned system
(ammonia–water); every other case — hydrocarbons, alcohol–water, acetone–water,
the extractive and azeotropic ternaries — is a neutral (non-electrolyte) system
and needs none of it. A Helmholtz EOS is either single-system (Tillner-Roth) or a
from-scratch reference backbone nothing planned requires (the cubic EOS already
meets the roughly 1 percent cross-simulator target, and the validation oracles
use cubic EOS themselves). NRTL, by contrast, is general infrastructure: it
improves every aqueous-organic column on the ladder and is the standard model for
the later extractive and azeotropic cases. So the honest engineering choice is to
add NRTL once (broad benefit) and reproduce the ammonia–water textbook chart from
reference data rather than build single-use thermodynamics.

**The two constructions in this notebook.** Route (a) computes the H–x–y diagram
from NRTL and steps off stages — teaching how a chart is built from thermo, and
showing the fray. Route (b) feeds reference enthalpy + VLE points through the
enthalpy-augmented `from_points` path and runs the same construction —
reproducing the textbook faithfully and cleanly separating the method from the
thermo model. Shown side by side, the gap between them is the lesson.
