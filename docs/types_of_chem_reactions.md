## Chemical Reaction Types

Here’s we list different types of chemical reactions using abstract species A, B, C, etc., including various forms like catalytic, reversible, and reactions with different stoichiometries.


### Simple Irreversible Reaction
This is the most basic type of reaction, where species \(A\) reacts with species \(B\) to form product \(C\).

**Reaction:**
$$
A + B \rightarrow C
$$

**Rate Law:**
$$
r = k[A][B]
$$

- \(k\) is the rate constant.
- The reaction rate depends on the concentrations of \(A\) and \(B\).

---

### Second-Order Reaction (Stoichiometric Coefficients)
Here, multiple moles of \(A\) react with one mole of \(B\).

**Reaction:**
$$
2A + B \rightarrow C
$$

**Rate Law:**
$$
r = k[A]^2[B]
$$

- The rate is proportional to the square of the concentration of \(A\) and the concentration of \(B\).

---

### Reversible Reaction
A reversible reaction allows products to convert back to reactants.

**Reaction:**
$$
A + B \rightleftharpoons C + D
$$

**Forward Rate Law:**
$$
r_{\text{forward}} = k_{\text{f}}[A][B]
$$

**Reverse Rate Law:**
$$
r_{\text{reverse}} = k_{\text{r}}[C][D]
$$

- \(k_{\text{f}}\) and \(k_{\text{r}}\) are the forward and reverse rate constants, respectively.
- At equilibrium, the forward and reverse rates are equal.

---

### Catalytic Reaction
In a catalytic reaction, a catalyst (denoted as \(X\)) speeds up the reaction without being consumed.

**Reaction:**
$$
A + B + X \rightarrow C + X
$$

**Rate Law:**
$$
r = k[A][B][X]
$$

- The catalyst \(X\) appears in the rate law but is not consumed in the overall reaction.
- Catalysts reduce the activation energy of the reaction but remain unchanged.

---

### Autocatalytic Reaction
In an autocatalytic reaction, one of the products of the reaction acts as a catalyst.

**Reaction:**
$$
A + B \rightarrow 2C
$$

**Rate Law:**
$$
r = k[A][B][C]
$$

- Here, the product \(C\) catalyzes the reaction, so the reaction rate increases as more product forms.

---

### Third-Order Reaction
A reaction involving three reactants in specific stoichiometric amounts.

**Reaction:**
$$
A + B + C \rightarrow D
$$

**Rate Law:**
$$
r = k[A][B][C]
$$

- This reaction is third-order overall since it involves three different reactants.

---

### Parallel Reactions
In parallel reactions, a single reactant can follow multiple pathways to form different products.

**Reaction:**
$$
A \rightarrow B \quad \text{and} \quad A \rightarrow C
$$

**Rate Laws:**
$$
r_B = k_1[A] \quad \text{and} \quad r_C = k_2[A]
$$

- \(A\) can convert into either \(B\) or \(C\), with different rate constants \(k_1\) and \(k_2\) for each pathway.

---

### Consecutive Reactions
In consecutive reactions, the product of one reaction becomes the reactant in another reaction.

**Reaction:**
$$
A \rightarrow B \rightarrow C
$$

**Rate Laws:**
$$
r_1 = k_1[A], \quad r_2 = k_2[B]
$$

- \(A\) first converts to \(B\), and then \(B\) converts to \(C\). The reaction rate of the second step depends on the concentration of \(B\).

---

### Inhibition Reaction
In an inhibition reaction, one species slows down or inhibits the reaction.

**Reaction:**
$$
A + B \rightarrow C \quad \text{(inhibited by \(D\))}
$$

**Rate Law:**
$$
r = \frac{k[A][B]}{1 + [D]}
$$

- The presence of \(D\) inhibits the reaction, reducing the reaction rate as the concentration of \(D\) increases.

---

### Enzymatic (Michaelis-Menten) Reaction
This is a specialized catalytic reaction where an enzyme \(E\) catalyzes the conversion of substrate \(S\) to product \(P\).

**Reaction:**
$$
S + E \rightleftharpoons SE \rightarrow P + E
$$

**Rate Law (Michaelis-Menten Approximation):**
$$
r = \frac{V_{\text{max}}[S]}{K_m + [S]}
$$

- $(V_{\text{max}})$ is the maximum reaction rate, and \(K_m\) is the Michaelis constant.
- The enzyme catalyzes the reaction without being consumed.

---

### Zero-Order Reaction
In a zero-order reaction, the rate of reaction is independent of the concentration of reactants.

**Reaction:**
$$
A \rightarrow B
$$

**Rate Law:**
$$
r = k
$$

- The reaction rate remains constant, regardless of the concentration of \(A\).

---

### Bimolecular Reaction (Different Stoichiometry)
This involves a reaction between two molecules, but with different stoichiometric coefficients.

**Reaction:**
$$
3A + 2B \rightarrow C
$$

**Rate Law:**
$$
r = k[A]^3[B]^2
$$

- The rate depends on the cubic concentration of \(A\) and the square concentration of \(B\).


## Composite and Parallel Reactions

It is possible to have a reaction where the rate follows the form \( r = k_1[A] + k_2[B] \), but this is **not a typical form** for a single-step elementary chemical reaction. This rate expression might describe a **composite** or **multi-step** reaction mechanism or even a **parallel reaction** where different reactants independently lead to different products.

### Breakdown of the Rate Expression

The rate law \( r = k_1[A] + k_2[B] \) suggests that two processes are happening in parallel, each involving a different reactant. This form could be observed in specific types of reactions:

#### **Parallel Reactions**

In parallel reactions, a single reactant (or set of reactants) can undergo different reactions, each leading to different products. For example:

$$
A \rightarrow P \quad \text{(with rate constant } k_1\text{)}
$$

$$
B \rightarrow Q \quad \text{(with rate constant } k_2\text{)}
$$

In this scenario, the two reactants \(A\) and \(B\) independently produce different products, and the overall reaction rate is the sum of the rates of the individual reactions. 

- The rate of disappearance of \(A\) follows \( r_1 = k_1[A] \).
- The rate of disappearance of \(B\) follows \( r_2 = k_2[B] \).

The **total rate** of the overall reaction would be:

$$
r = k_1[A] + k_2[B]
$$

This type of reaction is common in competitive or alternative pathways, where two different species react independently in parallel processes.

#### First-Order Reactions in Different Reactants

In cases where both \(A\) and \(B\) are involved in different first-order reactions, each of them affects the rate independently. For example, the system could involve two sequential reactions happening at the same time, where each reactant controls its part of the reaction rate.

#### Other types

- Non-Elementary Reactions
    - If the reaction is a **non-elementary reaction**, the rate law can take a form that is not directly related to the stoichiometric coefficients of a balanced chemical equation. In these cases, the rate law is often derived experimentally and might involve sums of different terms. For instance, a complex multi-step mechanism or a surface-mediated reaction could yield a rate expression of this form.
- Surface Reactions or Catalytic Reactions
    - In certain catalytic reactions, especially those occurring on surfaces (e.g., heterogeneous catalysis), the reaction rate might depend on the concentration of multiple species in a linear fashion. For example, if two reactants \(A\) and \(B\) interact with a surface in independent steps, the overall reaction rate might be described by a sum of two independent terms.


---

### Final Remarks

This list covers a variety of reaction types, including simple irreversible reactions, catalytic reactions, reversible reactions, and those with different stoichiometries. Each reaction type has its corresponding rate law that describes how the reaction proceeds over time based on the concentrations of the involved species. 

In this project we leverage parallel reactions and reactions of different types to implement mathematical equations with multiple terms.