# McCabe–Thiele — the equations as implemented

![McCabe–Thiele staircase for a benzene–toluene 95/5 split at 101.325 kPa: blue
equilibrium curve, gray y = x diagonal, red rectifying and green stripping
operating lines meeting at the purple q-line, the numbered blue staircase of
theoretical stages, and the dashed limiting line at minimum reflux.](img/mccabe-thiele.png)

*Above: the diagram this document builds, rendered straight from the library —
benzene–toluene, $x_D = 0.95$, $x_B = 0.05$, $z_F = 0.5$, saturated-liquid feed,
$R = 1.5\,R_{\min}$. It gives $R_{\min} = 1.163$ and $N = 12.22$ stages with the
feed on stage 6 (the pinned notebook values, §6). Regenerate with
`stages.mccabe_thiele` + `stages.plotting.plot_mccabe_thiele`.*

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
$(x_i,\,y^\ast_i,\,T_i)$; queries interpolate linearly (binary search +
piecewise-linear, `interp`). The endpoints $y^\ast(0)=0$, $y^\ast(1)=1$ are pinned
exactly (thermodynamic identity).

The textbook idealization, kept as a constructor and used as the analytic test
oracle:

$$y = \frac{\alpha x}{1 + (\alpha - 1)\,x} \tag{S\&H 7-13}$$

Point relative volatility from the curve:
$\alpha(x) = \dfrac{y/(1-y)}{x/(1-x)}$.

## 2. Material balances (`column/model.rs`)

$$F = D + B, \qquad F z_F = D x_D + B x_B
\;\Rightarrow\;
D = F\,\frac{z_F - x_B}{x_D - x_B} \tag{S\&H 7-2, 7-3}$$

## 3. Construction lines (`binary/mccabe_thiele.rs`)

With external reflux ratio $R = L/D$ and constant molal overflow (the
McCabe–Thiele assumption — relaxed at rung 2, Ponchon–Savarit):

$$\text{rectifying:}\quad y = \frac{R}{R+1}x + \frac{x_D}{R+1} \tag{S\&H 7-9}$$

$$\text{q-line:}\quad y = \frac{q}{q-1}x - \frac{z_F}{q-1},
\quad\text{vertical at } q = 1 \tag{S\&H 7-26}$$

where $q$ is the fraction of the feed joining the liquid ($\bar L = L + qF$,
$\bar V = V - (1-q)F$). The stripping line is taken in two-point form through
$(x_B, x_B)$ and the rectifying ∩ q-line intersection
(`operating_intersection`), which is equivalent to S&H eq. 7-12.

## 4. Stage stepping (`mccabe_thiele`, `total_reflux`)

From $(x_D, x_D)$, alternate: **horizontal** to the curve
($x_n$ such that $y^\ast(x_n) = y_n$, via inverse interpolation) — one theoretical
stage — then **vertical** to the active operating line
($y_{n+1} = y_\text{op}(x_n)$). The staircase switches from the rectifying to
the stripping line the first time $x_n$ passes the operating-line intersection
— that stage is the **optimal feed stage**. Stepping ends when $x_n \le x_B$;
the final stage is counted fractionally, linearly in $x$:

$$N = (n-1) + \frac{x_{n-1} - x_B}{x_{n-1} - x_n}$$

(Fenske's closed form measures the last fraction logarithmically in the
composition ratio — the two conventions differ by ≲ 0.1 stage, never a whole
stage; asserted in the tests.)

At **total reflux** both operating lines collapse onto $y = x$ and the same
stepping yields $N_{\min}$ — on a constant-α curve this reproduces Fenske:

$$N_{\min} = \frac{\ln\!\big[(x_D/(1-x_D))\,\cdot\,((1-x_B)/x_B)\big]}{\ln \alpha}$$

**Condenser kinds:** a total condenser is not an equilibrium stage (stepping
starts on the diagonal at $x_D$); a partial condenser is one (the first step of
the staircase *is* the condenser). A partial (kettle) reboiler is always the
last equilibrium stage.

**Murphree vapor efficiency** ($E_{MV} < 1$, S&H §7.4): the horizontal step
targets the pseudo-curve

$$y_\text{eff}(x) = y_\text{op}(x) + E_{MV}\,\big(y^\ast(x) - y_\text{op}(x)\big)$$

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

   $$R_\text{strip} = \frac{q + s(1-q) - s\,d}{d\,(s - 1)}$$

4. $R_{\min} = \max(R_\text{rect},\, R_\text{strip})$; the arg-max sample is
   the reported pinch, flagged `tangent` when it sits away from the feed point.

For a concave curve both sections give the same feed pinch (asserted in tests);
on constant-α curves with $q = 1$ the result matches Underwood's binary closed
form (S&H eq. 7-24):

$$R_{\min} = \frac{1}{\alpha - 1}\left[\frac{x_D}{z_F} -
\alpha\,\frac{1 - x_D}{1 - z_F}\right]$$

## 6. Validation summary

- Constant-α oracle tests: stepping vs Fenske, `rmin` vs Underwood
  (`binary/mccabe_thiele.rs` unit tests).
- Real-thermo integration: benzene–toluene (PR, α ≈ 2.3–2.5, boiling-point
  endpoints) and methanol–water (van Laar Λ₁₂ = 0.5853, Λ₂₁ = 0.3458 from vle
  Chapter IV / Orbey & Sandler) — `engine/tests/m1_mccabe_thiele.rs`.
- Notebook-pinned values (executed assertions): R_min = 1.163 and N = 12.22 at
  R = 1.5 R_min for the benzene–toluene 95/5 split; Table 4.6 bubble-pressure
  reproduction for methanol–water; ethanol–water tangent pinch at x ≈ 0.735.
