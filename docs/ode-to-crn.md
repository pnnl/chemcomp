# From polynomial ODEs to chemical reaction networks

ChemComp compiles a system of ordinary differential equations (ODEs) into an
abstract chemical reaction network (CRN): a set of reactions whose
**mass-action kinetics are exactly the given equations**. This page explains
the math behind the flow and the passes that implement it.

```
Modelica / bmodelica ──convert-bmodelica-to-ode──▶ ode.system
  ──ode-normalize──▶ ──ode-dual-rail──▶ ──convert-ode-to-acr──▶ acr.network
  ──acr-network-merge──▶ ──acr-network-to-crn──▶ acr.crn ──▶ crn-db-match / mlir-to-json
```

`chemcomp-opt model.mlir --ode-to-crn` runs everything after
`convert-bmodelica-to-ode`.

## 1. Mass-action kinetics

A reaction turns a *complex* of reactants into a complex of products, e.g.
`X + Y → 2Y`. Write a complex as a vector of counts, `y = (y_1, …, y_n)`, so
`X + Y` is `(1, 1)` and `2Y` is `(0, 2)`. Under mass action, the reaction
`y → y'` with rate constant `k` proceeds at speed `k · x^y`, where
`x^y = x_1^{y_1} ⋯ x_n^{y_n}`, and changes the amounts by `y' − y`. Summing over
the reactions of a network gives its ODEs:

$$\dot x = \sum_{y \to y'} k_{y \to y'}\, x^{y}\, (y' - y)$$

Every right-hand side is a polynomial in the concentrations. The question this
flow answers is the converse: **given a polynomial ODE system, find reactions
whose mass-action ODEs are that system.**

Example: Lotka–Volterra, `x' = a·x − b·x·y`, `y' = b·x·y − c·y`, is realized by

| Reaction | Speed | Contribution |
|---|---|---|
| `X → 2X` (rate `a`) | `a·x` | `+a·x` to `x'` |
| `X + Y → 2Y` (rate `b`) | `b·x·y` | `−b·x·y` to `x'`, `+b·x·y` to `y'` |
| `Y → ∅` (rate `c`) | `c·y` | `−c·y` to `y'` |

## 2. Which ODEs are possible: the Hungarian Lemma

Not every polynomial system has such a network. A reaction can only consume
`X` if `X` is one of its reactants, so its speed `k·x^y` then contains `x`.
Hence, in `x'`, **every negative term must contain `x`**. A negative term
without `x` is called a *negative cross-effect*.

**Hungarian Lemma** (Hárs & Tóth 1981). A polynomial ODE system is the
mass-action system of some reaction network if and only if it has no negative
cross-effects.

The condition is necessary by the argument above. It is sufficient because of
the *canonical realization*: each term `c·x^y` of `x_i'` becomes one reaction

- `c > 0`: `y → y + X_i`, with rate `c`;
- `c < 0`: `y → y − X_i`, with rate `|c|`. This is a valid reaction (no
  negative product count) exactly because `x^y` contains `x_i`.

The reaction `y → y ± X_i` changes only `x_i`, by `±1`, at speed `|c|·x^y`, so
it contributes exactly the term `c·x^y` to `x_i'` and nothing else. Summing
over all terms gives back the whole system.

Systems with cross-effects include the harmonic oscillator `x' = −y`, `y' = x`
(`−y` in `x'`) and Lorenz (`−x·z` in `y'`). Their solutions also leave the
non-negative orthant, where concentrations can't go.

`ode-check-kinetic` reports every cross-effect, e.g.
`der(x) has -y, which does not contain x`.

## 3. Dual-rail encoding

To compile systems with cross-effects, represent a state as the difference of
two non-negative concentrations, `x = x⁺ − x⁻` (`x_p`, `x_m` in the IR).
`ode-dual-rail`:

1. substitutes `x_j = x_j⁺ − x_j⁻` for each split state and expands, e.g.
   `x·z` becomes `x⁺z⁺ − x⁺z⁻ − x⁻z⁺ + x⁻z⁻`;
2. routes each resulting term of `x'` to one rail. Since
   `x' = (x⁺)' − (x⁻)'`, a term can go to `(x⁺)'` unchanged or to `(x⁻)'` with
   the opposite sign. It picks the rail where the term is harmless:
   - positive term: to `(x⁺)'` (production), or, if it contains `x⁻`, to
     `(x⁻)'` as a negative term (degradation of `x⁻`);
   - negative term: stays in `(x⁺)'` if it contains `x⁺` (degradation of
     `x⁺`), otherwise moves to `(x⁻)'` as a positive term (production);
3. adds `−η·x⁺x⁻` to both rails: the annihilation `X⁺ + X⁻ → ∅`.

After routing, every negative term of a rail contains that rail, so the
encoded system has no cross-effects and has a network.

**Exactness.** The annihilation term appears in both rails and cancels in
`(x⁺)' − (x⁻)'`, and the routed terms add up to the original right-hand side.
So `x⁺ − x⁻` follows the original equation exactly, for any `η > 0`. The rate
`η` only controls how large the rails get: without it both rails can grow
without bound while their difference stays correct. A large `η` keeps them
small but makes the system stiffer.

**Auto mode** (the default) splits only what is needed: the targets of the
cross-effects. Splitting `y` turns a term `y` elsewhere into `y⁺ − y⁻`, which
can create new cross-effects (in the oscillator, `x' = −y` becomes
`x' = −y⁺ + y⁻`), so the pass repeats until the system is kinetic. Split rails
never have cross-effects, so this takes at most one round per state.
`mode=all` splits every state; `states=…` adds states explicitly. Kinetic
systems (decay, Lotka–Volterra, …) are left alone.

Costs: two species per split state, a term of degree `d` in split states
expands to up to `2^d` terms, and one annihilation reaction per pair.

The oscillator `x' = −w·y`, `y' = w·x` becomes

```
x_p' = w·y_m − x_p·x_m     x_m' = w·y_p − x_p·x_m
y_p' = w·x_p − y_p·y_m     y_m' = w·x_m − y_p·y_m
```

i.e. four catalytic productions (`Y⁻ → Y⁻ + X⁺`, …) and two annihilations.

## 4. Composable rules

The canonical realization always works, but it gives one reaction per term,
and those reactions rarely look like chemistry: Lotka–Volterra's predation
becomes two half-reactions, `X + Y → Y` and `X + Y → X + 2Y`. ChemComp
therefore composes networks from **rules** first and uses the canonical
realization only for what no rule covers.

A rule is one reaction motif over rule-local species, defined in
`include/ChemComp/Conversion/ODEToACR/ODERules.td` and compiled into a table
(`chemcomp-tblgen -gen-ode-rules`). The set is fixed:

| Rule | Motif | Terms claimed |
|---|---|---|
| binding | A + B → C | 3 |
| dissociation | C → A + B | 3 |
| predation | X + Y → 2Y | 2 |
| cubic-autocatalysis | 2X + Y → 3X | 2 |
| conversion | A → B | 2 |
| annihilation | X + Y → ∅ | 2 |
| catalytic-production | Y → Y + X | 1 |
| autocatalysis | X → 2X | 1 |
| second-order-decay | 2X → X | 1 |
| degradation | X → ∅ | 1 |
| production | ∅ → X | 1 |

**Term signature.** By mass action, the motif `y → y'` with rate `k` adds
`k·(y'_s − y_s)·x^y` to `s'` for each species `s` whose count changes. For
predation, `X + Y → 2Y` adds `−k·x·y` to `x'` and `+k·x·y` to `y'`. Because the
signature is derived from the reaction, a rule can't disagree with what it
emits.

**Matching.** A rule matches when some injective assignment of its species to
states and some `k > 0` (a number times a product of parameters) make every
term of its signature appear in the system with exactly that coefficient. The
rule then *claims* those terms and emits its reaction. Rate constants are
compared symbolically: `−b1·x·y` and `+b2·x·y` are not one predation.

**Composition.** Rules are applied greedily by priority (larger motifs first)
until none matches, then:

- `canonical-fallback=true` (default): each unclaimed term becomes its
  canonical reaction. The flow therefore compiles every kinetic system;
- `canonical-fallback=false`: unclaimed terms are an error listing them, for
  when only rule-made reactions are acceptable.

Since every emitted reaction contributes exactly the terms it claimed, and
every term is claimed once, the network's ODEs are the system's ODEs.

Each reaction records its `rule`, its rate expression (`rate_name`, e.g.
`"b"`, `"0.5*k"`), its numeric `rate` when all parameters have values, and the
terms it accounts for (`from`):

```mlir
acr.rxn @x + @y -> 2 @y rate 0.4 {from = "der(x): -b*x*y; der(y): b*x*y", rate_name = "b", rule = "predation"}
```

## 5. Merging reactions

Different networks can have the same ODEs. What matters for a reactant
complex `y` is only its *net reaction vector* `Σ k (y' − y)` over its reactions
(Buxton et al. 2025, eq. 2–3). `acr-network-merge` uses this:

- identical reactions are summed (their rates add);
- two reactions from the same complex at the same rate are one reaction with
  the summed change, if no product count becomes negative. For example,
  `X + Y → Y` and `X + Y → X + 2Y`, both at rate `b`, become `X + Y → 2Y`;
- reactions whose products equal their reactants are dropped.

It mostly tidies up after the canonical fallback.

Finding the "best" network among all realizations is a research topic of its
own. For example, Buxton et al. (2025) decide whether a weakly reversible,
deficiency-zero realization exists, which guarantees stable dynamics. That is
not implemented here.

## 6. Passes and tools

| Pass / tool | What it does |
|---|---|
| `convert-bmodelica-to-ode` | Base Modelica model → `ode.system` (needs MARCO) |
| `ode-normalize` | sum like terms, canonical order |
| `ode-check-kinetic` | report negative cross-effects (`require` makes them errors) |
| `ode-dual-rail` | `mode=auto\|all`, `states=…`, `annihilation-rate=η` |
| `convert-ode-to-acr` | rules + canonical fallback (`canonical-fallback=false`) |
| `acr-network-merge` | merge reactions without changing the ODEs |
| `acr-network-to-crn` | `acr.network` → `acr.crn` / `acr.crn_generic` |
| `ode-to-crn` | the pipeline; options `dual-rail=auto\|all\|none`, `annihilation-rate`, `canonical-fallback`, `merge` |
| `chemcomp-translate -mlir-to-json` | network → emulator JSON (with catalysts, initial amounts, `species_map`) |
| `chemcomp-translate -ode-to-json` | `ode.system` → reference ODE in JSON |

```sh
chemcomp-opt model.mlir --ode-to-crn                                  # network as acr.crn
chemcomp-opt model.mlir --ode-to-crn=canonical-fallback=false         # rules only
chemcomp-opt model.mlir --ode-normalize --ode-dual-rail --convert-ode-to-acr \
    --acr-network-merge                                               # readable acr.network
chemcomp-opt model.mlir --convert-bmodelica-to-ode --ode-to-crn \
  | chemcomp-translate -mlir-to-json > network.json                   # from Modelica
```

## 7. Checking a network numerically

`tools/chemcomp-reservoir-emulator/scripts/check_realization.py` solves the
reference ODE, simulates the network with the emulator's mass-action
dynamics, maps rails back (`x = x_p − x_m`, from `species_map`), and reports
the mean squared error per state. `pixi run check-odes` runs it on the
examples; their lit tests run it too when the emulator is installed.

| Example | States | Terms | Degree | Kinetic | Network | Worst MSE |
|---|---|---|---|---|---|---|
| decay | 1 | 1 | 1 | yes | X → ∅ | 0 |
| birth–death | 1 | 2 | 0–1 | yes | ∅ → X, X → ∅ | 0 |
| dimerization `x' = −k·x²` | 1 | 1 | 2 | yes | 2X → X | 1e-32 |
| reversible A ⇌ B | 2 | 4 | 1 | yes | A → B, B → A | 0 |
| Lotka–Volterra | 2 | 4 | 2 | yes | X → 2X, X + Y → 2Y, Y → ∅ | 0 |
| SIR | 3 | 4 | 2 | yes | S + I → 2I, I → R | 6e-34 |
| Michaelis–Menten | 4 | 9 | 2 | yes | E + S → C, C → E + S, C → E + P | 0 |
| Brusselator | 2 | 6 | 0–3 | yes | ∅ → X, X → Y, 2X + Y → 3X, X → ∅ | 2e-19 |
| harmonic oscillator | 2 | 2 | 1 | no | dual-rail: 4 catalytic productions, 2 annihilations | 1e-20 |
| Lorenz (t ≤ 1) | 3 | 7 | 2 | no | dual-rail: 13 rule reactions + 8 canonical | 7e-17 |

The errors are integrator round-off; the networks are exact. Lorenz uses the
canonical fallback for its bimolecular catalytic terms (e.g. `x⁺·z⁻` producing
`y⁺`), which no rule covers; with `canonical-fallback=false` it fails and
names them. The sources are in `test/Conversion/ODEToACR/examples/`.

## 8. Limitations

- Only explicit polynomial ODEs: one `der(x) = polynomial` per state.
  Non-polynomial terms (`sin`, division by a state) are rejected; recasting
  them into polynomial form (quadratization) is future work.
- Parameters are rate constants and must be non-negative; signs belong in the
  equations.
- Dual-rail networks contain `X⁺ + X⁻ → ∅`, which a chemical database is
  unlikely to provide; `crn-db-match` works best on kinetic models.
- The greedy rule order is fixed. It doesn't search for the smallest network,
  and the rule set can't be extended at runtime.

## References

- V. Hárs, J. Tóth. *On the inverse problem of reaction kinetics.* In
  Qualitative Theory of Differential Equations (Szeged, 1979), Colloquia
  Mathematica Societatis János Bolyai 30, 1981, 363–379.
- P. Érdi, J. Tóth. *Mathematical Models of Chemical Reactions.* Manchester
  University Press, 1989.
- N. Buxton, G. Craciun, A. Deshpande, C. Pantea. *Weakly reversible
  deficiency zero realizations of reaction networks.* Journal of Nonlinear
  Science 35:104, 2025.
- K. Oishi, E. Klavins. *Biomolecular implementation of linear I/O systems.*
  IET Systems Biology 5(4), 2011 (dual-rail representation of signed signals).
