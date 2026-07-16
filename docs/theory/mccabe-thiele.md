# McCabe–Thiele — the equations as implemented

![McCabe–Thiele staircase for a benzene–toluene 95/5 split at 101.325 kPa: blue
equilibrium curve, gray y = x diagonal, red rectifying and green stripping
operating lines meeting at the purple q-line, the numbered blue staircase of
theoretical stages, and the dashed limiting line at minimum reflux.](img/mccabe-thiele.png)

**The diagram above** is rendered straight from the library — benzene–toluene,
$x_D = 0.95$, $x_B = 0.05$, $z_F = 0.5$, saturated-liquid feed, $R = 1.5 R_{\min}$.
It gives $R_{\min} = 1.163$ and $N = 12.22$ stages with the feed on stage 6 (the
pinned notebook values, §7). Regenerate with `stages.mccabe_thiele` +
`stages.plotting.plot_mccabe_thiele`.

*Modules: `engine/src/binary/equilibrium.rs`, `engine/src/binary/mccabe_thiele.rs`,
`engine/src/column/model.rs`. Notebook: `notebooks/01-mccabe-thiele.ipynb`.*

Primary sources: McCabe, W. L.; Thiele, E. W. Graphical Design of Fractionating
Columns. *Ind. Eng. Chem.* **1925**, 17 (6), 605–611. Implementation-level
treatment and equation numbers: Seader, J. D.; Henley, E. J.; Roper, D. K.
*Separation Process Principles* (S&H), Ch. 7.

Conventions: compositions are mole fractions of the **light** (more volatile)
component; stages are numbered **top-down** (1 = top; the partial reboiler is
the last stage). Units: K, kPa absolute, kmol/h.

## 1. The equilibrium curve (`binary/equilibrium.rs`)

At fixed column pressure $P$, each liquid composition $x$ has a bubble point
$T(x)$ with incipient vapor $y^\ast(x)$. `EquilibriumCurve::from_thermo` sweeps
`ThermoSystem::bubble_temperature` (vle-thermo) over an even $x$ grid and stores
$(x_i, y^\ast_i, T_i)$; queries interpolate linearly (binary search +
piecewise-linear, `interp`). The endpoints $y^\ast(0)=0$, $y^\ast(1)=1$ are pinned
exactly (thermodynamic identity).

The textbook idealization (S&H eq. 7-13), kept as a constructor and used as the
analytic test oracle:

$$y = \frac{\alpha x}{1 + (\alpha - 1) x}$$

Point relative volatility from the curve:
$\alpha(x) = \dfrac{y/(1-y)}{x/(1-x)}$.

## 2. Material balances (`column/model.rs`)

Overall and light-component balances (S&H eqs. 7-2, 7-3):

$$F = D + B, \qquad F z_F = D x_D + B x_B
\quad\Rightarrow\quad
D = F \frac{z_F - x_B}{x_D - x_B}$$

## 3. Construction lines (`binary/mccabe_thiele.rs`)

With external reflux ratio $R = L/D$ and constant molal overflow (the
McCabe–Thiele assumption — relaxed at rung 2, Ponchon–Savarit), the rectifying
line (S&H eq. 7-9) and the q-line (S&H eq. 7-26) are:

$$\text{rectifying:}\quad y = \frac{R}{R+1}x + \frac{x_D}{R+1}$$

$$\text{q-line:}\quad y = \frac{q}{q-1}x - \frac{z_F}{q-1},
\quad\text{vertical at } q = 1$$

where $q$ is the fraction of the feed joining the liquid ($\bar L = L + qF$,
$\bar V = V - (1-q)F$). The stripping line is taken in two-point form through
$(x_B, x_B)$ and the rectifying ∩ q-line intersection
(`operating_intersection`), which is equivalent to S&H eq. 7-12.

## 4. Stage stepping (`mccabe_thiele`, `total_reflux`)

From $(x_D, x_D)$, alternate: **horizontal** to the curve
($x_n$ such that $y^\ast(x_n) = y_n$, via inverse interpolation) — one theoretical
stage — then **vertical** to the active operating line, giving
$y_{n+1} = y_\text{op}(x_n)$. The staircase switches from the rectifying to
the stripping line the first time $x_n$ passes the operating-line intersection
— that stage is the **optimal feed stage**. Stepping ends when $x_n \le x_B$;
the final stage is counted fractionally, linearly in $x$:

$$N = (n-1) + \frac{x_{n-1} - x_B}{x_{n-1} - x_n}$$

(Fenske's closed form measures the last fraction logarithmically in the
composition ratio — the two conventions differ by ≲ 0.1 stage, never a whole
stage; asserted in the tests.)

At **total reflux** both operating lines collapse onto $y = x$ and the same
stepping yields $N_{\min}$ — on a constant-α curve this reproduces Fenske:

$$N_{\min} = \frac{\ln\big[(x_D/(1-x_D)) \cdot ((1-x_B)/x_B)\big]}{\ln \alpha}$$

**Condenser kinds:** a total condenser is not an equilibrium stage (stepping
starts on the diagonal at $x_D$); a partial condenser is one (the first step of
the staircase *is* the condenser). A partial (kettle) reboiler is always the
last equilibrium stage.

**Murphree vapor efficiency** ($E_{MV} < 1$, S&H §7.4): the horizontal step
targets the pseudo-curve

$$y_\text{eff}(x) = y_\text{op}(x) + E_{MV} \big(y^\ast(x) - y_\text{op}(x)\big)$$

inverted by bisection ($y_\text{eff}$ is monotone). The pseudo-curve is applied
to every stage including the reboiler — slightly conservative (a real reboiler
is a true equilibrium stage); documented in the module. $R_{\min}$ always uses
the true curve (the classical construction).

## 5. Minimum reflux and pinch detection (`rmin`)

The pinch is found geometrically over the sampled curve — **not** assumed to
sit at the feed point, so tangent pinches are handled on both sections:

1. **Feed point**: the q-line ∩ curve intersection $(x_q^\ast, y_q^\ast)$
   (`q_line_curve_intersection`; special-cased vertical/horizontal, otherwise a
   sign-change scan outward from $z_F$, on the side determined by the q-line
   slope).
2. **Rectifying side**: the operating line anchored at $(x_D, x_D)$ must stay
   below the curve on $[x_q^\ast, x_D]$. Limiting slope
   $m = \max_e \dfrac{x_D - y_e}{x_D - x_e}$ over curve samples (the feed point
   is a candidate), then $R_\text{rect} = m/(1-m)$. $m \ge 1$ ⇒ the spec is
   **infeasible** at any reflux (e.g. $x_D$ beyond an azeotrope) —
   `StagesError::Infeasible`.
3. **Stripping side**: the line anchored at $(x_B, x_B)$ must stay below the
   curve on $[x_B, x_q^\ast]$; limiting slope
   $s = \min_e \dfrac{y_e - x_B}{x_e - x_B}$. With $d = D/F$ from §2, the
   feed-section balances convert $s = \bar L/\bar V$ to an equivalent reflux:

   $$R_\text{strip} = \frac{q + s(1-q) - s d}{d (s - 1)}$$

4. $R_{\min} = \max(R_\text{rect}, R_\text{strip})$; the arg-max sample is
   the reported pinch, flagged `tangent` when it sits away from the feed point.

For a concave curve both sections give the same feed pinch (asserted in tests);
on constant-α curves with $q = 1$ the result matches Underwood's binary closed
form (S&H eq. 7-24):

$$R_{\min} = \frac{1}{\alpha - 1}\left[\frac{x_D}{z_F} -
\alpha \frac{1 - x_D}{1 - z_F}\right]$$

## 6. Assumptions and limitations — what the diagram takes for granted

McCabe–Thiele is rung 1 of the ladder because it buys its clarity with
assumptions. As implemented here they are:

1. **Constant molal overflow (CMO)** — the load-bearing one, unpacked below.
2. **Constant column pressure** — a single $P$ for every stage; the equilibrium
   curve is sampled once at that pressure, with no tray-to-tray pressure drop.
3. **Theoretical (equilibrium) stages** — vapor and liquid leaving a stage are
   in exact equilibrium; real trays enter only through a uniform Murphree
   correction (§4), not stage-by-stage efficiencies.
4. **Binary mixtures only** — compositions live on one axis; a third component
   has nowhere to go (rungs 3–4 lift this).
5. **The feed's thermal state is one number, $q$** — heat losses, interstage
   heaters/coolers and multiple feeds are out of scope at this rung
   (Ponchon–Savarit handles side duties naturally).
6. **No energy answers.** The diagram never touches an enthalpy, so it cannot
   produce $Q_C$ or $Q_R$. Sizing a reboiler from rung 1 means bolting on a
   separate $Q \approx V \lambda$ estimate — which inherits every CMO error
   below.

### 6.1 What "constant molal overflow" actually says

CMO is the claim that, within each column section, the **molar** flows do not
change from stage to stage:

$$L_0 = L_1 = \dots = L \quad\text{and}\quad V_1 = V_2 = \dots = V
\quad\text{(rectifying; barred symbols for stripping)}$$

which is exactly what makes the operating lines of §3 **straight**: with $L$
and $V$ constant, the component balance around the top of the column,
$y_{n+1} = (L/V) x_n + (D/V) x_D$, has a constant slope $L/V$.

CMO is an **energy statement in disguise**. Write the energy balance around
stage $n$ of an adiabatic column section (S&H §7.2):

$$V_{n+1} H_{V,n+1} + L_{n-1} h_{L,n-1} = V_n H_{V,n} + L_n h_{L,n}$$

Three idealizations make the enthalpies drop out of it:

1. **equal molar latent heats** — both components carry the same
   $\lambda$ per kmol of phase change;
2. **negligible sensible heat** — the section's temperature span is small
   enough that $c_p \Delta T \ll \lambda$;
3. **no heat of mixing** — $h^E \approx 0$, so a stream's enthalpy is just the
   mole-weighted sum of its components'.

Under those three, $H_{V,n} \approx h_L(T_n) + \lambda$ for every stage, all
the $h_L$ terms cancel, and the balance collapses to $V_{n+1} = V_n$; the mass
balance then forces $L_n = L_{n-1}$. The physical picture is a **one-for-one
molar trade**: condensing one kmol of vapor on a stage releases $\lambda$,
which is exactly the heat needed to boil up one kmol of liquid — so whatever
vapor arrives, the same molar flow leaves. This is also why distillation
bookkeeping is done in **moles**, not mass: by Trouton's rule
($\lambda \approx 88 \cdot T_b$ kJ/kmol), chemically similar, close-boiling
pairs really do have nearly equal *molar* latent heats, while their *mass*
latent heats can differ wildly. CMO in mass units would fail even for
benzene–toluene.

### 6.2 When it fails, and how badly

Each idealization above maps to a failure mode: unequal $\lambda$ (wide-boiling
or chemically dissimilar pairs), a large temperature span (sensible heat), and
associating/aqueous systems (heat of mixing). The energy-exact construction on
the enthalpy–composition diagram — [Ponchon–Savarit](ponchon-savarit.md), rung
2 — shows what the flows really do:

![Internal molar flows per mole of feed versus liquid composition: for
benzene–toluene the energy-exact liquid and vapor profiles sit within a few
percent of the flat CMO lines; for ammonia–water they run 25–40% below
them.](img/cmo-flow-profiles.png)

For benzene–toluene ($\lambda$ ratio 1.10) the true profiles hug the dashed CMO
levels to within a few percent — which is why the M2 consistency check agrees
across methods to within a stage. For ammonia–water
($\lambda$ ratio 1.74 plus a large exothermic heat of mixing) the internal
reflux runs 40% below what CMO assumes. And because the passing-stream ratio
$L/V$ is the **slope of the operating relation**, a varying $L/V$ bends the
straight operating line into a curve:

![x–y diagram for acetone–water: the true, energy-balanced operating curve
bends toward the equilibrium curve and pinches where the straight CMO
operating lines still show clearance; inset numbers show McCabe–Thiele
predicting 9.9 stages where the energy-exact answer is
14.0.](img/cmo-operating-lines.png)

The bite comes **near the pinch**, where the clearance between operating and
equilibrium curves is the whole game. For acetone–water ($\lambda$ ratio 1.31,
zero heat of mixing in the van Laar model — the failure is latent heats alone)
at $x_D = 0.90$, $x_B = 0.05$, $z_F = 0.30$:

- McCabe–Thiele reports $R_{\min} = 0.225$; the energy-exact minimum is
  $0.248$ — **CMO underestimates minimum reflux by 10%**, so a column designed
  at the classic heuristic $R = 1.1 R_{\min}^{CMO}$ **cannot meet spec at
  all** (the real construction pinches).
- At $R = 1.2 R_{\min}^{CMO}$ McCabe–Thiele counts 9.9 stages; the energy
  balance needs 14.0 — a **41% under-design**.

Both numbers are pinned as executed assertions in
[`notebooks/02-ponchon-savarit.ipynb`](../../notebooks/02-ponchon-savarit.ipynb)
(worked example 3). The rule of thumb: trust CMO for close-boiling, chemically
similar, near-ideal pairs at comfortable reflux; the moment the latent heats
diverge, the boiling range widens, the mixture associates, or the design
presses toward $R_{\min}$, move up a rung to the energy-exact method — that is
what it is for.

## 7. Validation summary

- Constant-α oracle tests: stepping vs Fenske, `rmin` vs Underwood
  (`binary/mccabe_thiele.rs` unit tests).
- Real-thermo integration: benzene–toluene (PR, α ≈ 2.3–2.5, boiling-point
  endpoints) and methanol–water (van Laar Λ₁₂ = 0.5853, Λ₂₁ = 0.3458 from vle
  Chapter IV / Orbey & Sandler) — `engine/tests/m1_mccabe_thiele.rs`.
- Notebook-pinned values (executed assertions): R_min = 1.163 and N = 12.22 at
  R = 1.5 R_min for the benzene–toluene 95/5 split; Table 4.6 bubble-pressure
  reproduction for methanol–water; ethanol–water tangent pinch at x ≈ 0.735.
