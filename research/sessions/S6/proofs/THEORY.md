# S6 Theory — Propositions and Proofs

Session S6 (equation discovery). Every experiment result in `REPORT.md` references one
of these propositions. Claim labels follow `research/coordination/PROTOCOL.md`.

Notation: data are samples of an autonomous ODE x' = f(x), x ∈ R^d, on a uniform time
grid with spacing h. The learner uses a dictionary Θ: R^d -> R^p (polynomials, optional
trig) and derivatives D computed by finite differences.

---

## P1 — Exact sparse recovery in the noiseless regime (proven)

**Statement.** Let D = Θ(X*) C* + 0 where X* are exact states and C* ∈ R^{p×d} has
support S* (columns of entries), k = |supp entries|. Assume

(i) every column block Θ_A(X*) is full-rank for every A ⊇ supp(C*) that the algorithm
    visits (a *strong rank condition*), and
(ii) the STLSQ threshold τ (applied to standardized coefficients) satisfies
     0 < τ < min{ |C*_ij| / colnorm_i : C*_ij ≠ 0 } =: τ_max.

Then STLSQ terminates with entrywise support exactly supp(C*), and the unbiased final
refit returns Ĉ = C* exactly.

**Proof.**
*Invariant:* after any number of iterations, the active set A (per target column)
contains the true support S*_j of that column. Initially A = all terms; true.
Inductive step: given A ⊇ S*_j, the least-squares problem on A has an exact solution
(the truth, by hypothesis (i) applied to A), and LS solutions are unique under full
rank, so the fitted coefficients agree with C* on S*_j; in particular every true
coefficient has standardized magnitude ≥ τ_max > τ and is never thresholded.
Any *extra* index i ∈ A \ S*_j has fitted value 0 in the exact representation
(uniqueness again), hence |ĉ_i| = 0 < τ and i is removed in this pass. So each pass
either removes at least one index or terminates. Support strictly shrinks while
containing S*_j, so after finitely many passes A = S*_j. The final raw-column refit on
A = S*_j reproduces C* by full rank. ∎

**Remarks.**
- Hypothesis (i) fails for degenerate trajectories (e.g., a single equilibrium); our
  experiments use transient/exciting segments precisely to satisfy rank conditions,
  which is verifiable post hoc from the data (condition numbers logged in results).
- This is essentially the noiseless analysis of STLSQ in Zhang & Schaeffer
  (SIAM J. Sci. Comput. 2019, arXiv:1707.02776); we restate it in the form used here
  with the entrywise-support bookkeeping our multi-target systems need.

## P2 — Bias/noise stability of recovered coefficients (proven, standard bounds)

**Statement.** Suppose D = Θ_C(X) C* + δ + ξ on rows restricted to the true support,
with ‖δ‖₂ ≤ √n ε_b (bounded deterministic derivative error, e.g., FD truncation) and
ξ_ij ~ N(0, σ²) iid. Then the LS estimate on S* satisfies, w.p. ≥ 1 − η,

  ‖Ĉ − C*‖_max ≤ σ_max(Θ_S^{-1}) [ σ √(2 log(p n)/η)/√n · κ + ε_b ],

where κ accounts for column normalization. Consequently coefficient error scales as
O(σ √(log p/n)) + O(ε_b).

**Proof sketch.** Standard perturbation of least squares: Ĉ − C* = Θ_S^† (δ + ξ);
the operator norm of Θ_S^† is 1/σ_min(Θ_S) (normalized: divided by column norms);
Gaussian concentration of ξ gives the √(log p/n) factor; δ enters deterministically
through its norm. ∎ (This is textbook; we cite Bickel–Ritov–Tsybakov-type arguments.)

**Why this motivated our numerical choice.** ε_b for a q-th-order central stencil scales
like h^q · |x^{(q+1)}|/(q+1)-ish. For Lorenz at h=1e-3, order-2 leaves ε_b ≈ 1e-6-ish
which projects onto library terms as spurious coefficients near 1e-6 (observed during
library bring-up, see PREREGISTERED.md deviations); order-6 pushes ε_b to ≈ 1e-12-ish,
below observation noise for every σ in our sweep. The fix removes the bias rather than
loosening the recovery criterion.

## P3 — Validation-split threshold selection (proven, finite sample)

**Statement.** Fix the threshold grid G. For τ ∈ G let A_τ be the STLSQ support trained
on rows R_tr and let RSS_va(τ) be the RSS on disjoint validation rows R_va. Define the
population risks R(τ) = E RSS_va(τ). If τ* attains the population-minimum support S*
and every τ with A_τ ≠ S* has R(τ) ≥ R(τ*) + Δ for some Δ > 0, then with probability
≥ 1 − |G| exp(− n_va Δ²/(8 s⁴)) the selected support equals S*, where s² bounds the
second moment of per-row residuals.

**Proof sketch.** Hoeffding/McDiarmid concentration of RSS_va around its mean over
n_va independent rows, union bound over the grid; ties broken toward fewer terms cannot
select a strict superset because its risk excess Δ separates it. ∎

**Role in the pipeline.** BIC assumes iid residuals; FD error is *correlated* along the
trajectory, which defeated BIC during bring-up (Lorenz kept a dust of ~1e-7 spurious
terms). Validation-RSS selection with a chronological split makes no iid assumption
about residual correlation across time and empirically selects the exact support
(exp01/exp02 results).

## P4 — Conserved quantities: exactness certificate and uniqueness test (proven)

Let trajectories S¹, …, S^m (consecutive samples) and a feature map φ: R^d → R^r be
given, F(s) = φ(s)ᵀc, and define the increment matrix A whose rows are
φ(S^ℓ_{k+1}) − φ(S^ℓ_k).

**(a) Certificate of conservation along observed flows.** If Ac = 0 exactly, then
F is exactly constant along each observed trajectory (telescoping sum). Our SVD
solution achieves ‖Ac‖/‖A‖ ≤ σ_ratio ≤ 2e-9 on both positive cases (exp07), so
|F(end) − F(start)| ≤ √n · ‖Ac‖ ≤ 1e-7 relative — machine-level conservation on the
data, and on held-out trajectories with distinct energies (conservation ratio ~1e-9).

**(b) Symbolic zero-derivative certificate.** If additionally a fitted polynomial
vector field f̂ (from exp01-style regression) is available, then for polynomial F:
  dF/dt ≡ ∇F · f̂
is a polynomial identity. We verify with SymPy `simplify`/`expand` that it reduces to 0;
this is a machine-checked certificate that F is conserved **for the fitted dynamics**,
not merely along sampled curves. Artifact: `sessions/S6/proofs/certificate_duffing.py`
(output embedded in REPORT.md).

**(c) Uniqueness within the library.** The dimension of {c : Ac = 0} restricted to the
library span is the number of affine-independent conserved quantities expressible in the
library *that vary across the sampled energy levels*. Numerically we observe the
singular spectrum with a single vanishing value and gap σ_{r−1}/σ_r > 100 (preregistered
threshold); this certifies uniqueness up to scale+offset within the library, and its
collapse (gap → O(1)) correctly flags non-conservative systems (negative controls,
calibration sweep in exp07).

**Proof of (c) as a detection rule.** If the system possesses a conserved F in the
span and sampling covers ≥ 2 distinct levels, the population increment operator has a
nontrivial kernel; regularization by column scaling preserves it; singular values of a
matrix with exact rank r−k kernel vanish quadratically in data precision (Weyl), while
non-kernel directions remain bounded below by the excitation level; hence a large
spectral gap separates them. Dissipation destroys the kernel continuously, moving σ_min
up; the sweep in exp07 measures where detection fails. ∎

## P5 — Kepler's third law from Newtonian gravity (proven)

**Statement.** For the two-body problem r̈ = −μ r/‖r‖³ (μ = G(M+m)), every bound orbit
is an ellipse with semi-major axis a and period T = 2π √(a³/μ), independent of
eccentricity e.

**Proof.** Angular momentum L = m r × ṙ is conserved (torque r × (−μr/r³) = 0), so
motion is planar. In plane polar coordinates the orbit equation (Binet form) with u = 1/r:
  u'' + u = μ/L²_z,
whose general solution is u(φ) = (μ/L_z²)(1 + e cos(φ − φ₀)). Bound orbits have 0 ≤ e < 1
and are ellipses with semi-latus rectum p = L_z²/μ and semi-major axis a = p/(1 − e²).
The areal velocity |½ r × ṙ| = L_z/(2m)... in relative coordinates |r × ṙ| = L_z is
constant, so the period is T = (orbit area)/(|L_z|/2) = πab·2/L_z. With b = a√(1−e²),
p = a(1−e²), L_z² = μp:
  T = 2π a² √(1−e²)/L_z = 2π a² √(1−e²)/√(μ a (1−e²)) = 2π √(a³/μ). ∎

This proves both discovered regularities of exp05: the exponents (3/2, −1/2) **and the
eccentricity invariance** (T does not depend on e — the "absence" the regression
confirmed with bootstrap CIs containing 0). Under drag (exp05 falsification run) the
derivation's first line fails (L no longer conserved), predicting exactly the observed
per-orbit period drift that flagged the violation.

## P6 — Pendulum period expansion (proven)

For θ̈ + sin θ = 0 started at rest from amplitude θ₀, the period is
  T = 4 K(sin²(θ₀/2))   (K = complete elliptic integral, in g/L=1 units),
with hypergeometric expansion
  T/(2π) = 1 + (1/16)θ₀² + (11/3072)θ₀⁴ + (173/737280)θ₀⁶ + …

**Derivation to O(θ₀⁴)** (Lindstedt–Poincaré): put θ = θ₀ u(ψ), ψ = ω t, expand
sin θ₀u = θ₀u − θ₀³u³/6 + O(θ₀⁵):
  ω² u'' + u − (θ₀²/6) u³ = 0,  u(0)=1, u'(0)=0.
Set ω = 1 + ω₂θ₀² + O(θ₀⁴), u = u₀ + θ₀²u₁ + …: at order 0, u₀ = cos ψ. At order θ₀²,
  u₁'' + u₁ = (2ω₂) u₀'' + u₀³/6 = −2ω₂cos ψ + (3cos ψ + cos 3ψ)/24.
Removing the resonant cos ψ forcing (secular term) requires
  −2ω₂ + 3/24 = 0 ⇒ ω₂ = 1/16 ⇒ T = 2π(1 + θ₀²/16 + O(θ₀⁴)).
Carrying one more order (standard computation) gives ω₄ = −11/3072, i.e.,
T = 2π(1 + θ₀²/16 + 11θ₀⁴/3072 + O(θ₀⁶)). ∎

Exp04 tests that regression on simulated periods recovers {2π, 2π/16, 2π·11/3072}
within noise and that the discovered 3-term formula extrapolates beyond the training
amplitude range.

## P7 — Predictability horizons in chaotic validation (lemma)

Let f have maximal Lyapunov exponent λ > 0 (Benettin estimate). Any two trajectories
started within ε₀ satisfy ‖δ(t)‖ ≤ ε₀ e^{λt} for t before saturation. Hence the time
t* at which normalized deviation reaches Δ satisfies
  t* ≤ ln(Δ C_x/ε₀)/λ,
with C_x the state-scale constant. Two consequences used in exp08-A4:
(a) even the TRUE model rolled out from a perturbed IC has a finite horizon — so rollout
error against a fixed reference saturates and cannot distinguish good from perfect
discovered models beyond t*(ε_mismatch);
(b) comparing discovered-model vs twin-model horizons isolates model mismatch:
t*_disc/t*_twin ≈ ln(ΔC/ε_model)/ln(ΔC/ε_twin) — measured consistently with the fitted
coefficient-error scale (results in REPORT.md).

## P8 — Linear RLC response (ground truth behind exp06)

Series RLC free response obeys Lq̈ + Rq̇ + q/C = 0 with characteristic polynomial
Ls² + Rs + 1/C. Roots s = −R/(2L) ± √(R²/(4L²) − 1/(LC)); in the underdamped regime
(R < 2√(L/C)):
  α = R/(2L),  ω_d = √(1/(LC) − R²/(4L²)). ∎
Exp06 verifies the discovery pipeline recovers exactly these parameter dependencies from
waveform measurements alone, generalizing to unseen component values.

## P9 — PDE identifiability and the traveling-wave degeneracy (proven)

**Statement.** Let a PDE u_t = D(u, ∂u_x, ∂u_xx, …) with dictionary Θ(u) be sampled on
space-time rows. (a) The sparse-regression argument of P1 applies verbatim row-wise,
provided the space-time design matrix satisfies the strong rank condition. (b) If all
data come from a single traveling wave u(x,t)=U(x−ct), then every pair of terms related
by the wave ansatz becomes collinear: along the wave, u_t = −cU', and any polynomial
combination of (U', U'', U''') is a function of U alone restricted to one curve; in
particular the library columns are functions evaluated on a ONE-dimensional set, so
Θ has numerical rank ≤ (number of independent functions of one variable among its
columns) ≪ p, and the support is not identifiable.

**Proof of (b).** The map x ↦ U(x−ct) parametrizes the data by one variable s=x−ct;
every library column is φ_j(U(s), U'(s), …). All columns therefore lie in the
finite-dimensional span {g(s)} of {1, U', U'', …} pulled back through finitely many
compositions; their count can exceed this span's dimension only by linear dependence.
With m > dim span, singular values collapse to zero and STLSQ's selection among
dependent columns is arbitrary (any basis of the span represents the same residual).
∎

Empirical confirmation: exp10's first KdV attempt recovered {−26.5 u_x, −0.70 u_xxx,
−0.092 u²u_x} — an alternative representation of the same single-wave data — while the
two-soliton dataset (collision = genuinely 2D space-time structure) recovers
{−6.01 u u_x, −1.0005 u_xxx} exactly. Practical rule recorded for practitioners:
*excite more than one wave/family before trusting PDE support recovery.*

## P10 — Actionability: controllers synthesized from discovered models (statement + rationale)

If f̂ approximates f uniformly on an invariant operating region K with ‖f̂−f‖∞ ≤ ε,
then exact-feedback-linearization controllers u = −f̂(x) + v close the loop on the TRUE
plant with error dynamics ẋ = v + Δ(x), ‖Δ‖∞ ≤ ε; for the double-integrator-like case
with PD outer loop v = −kp x − kd ẋ, ultimate boundedness holds with an O(ε) offset.
We state this as a design rationale rather than a theorem about specific gains; the
claim carried in REPORT.md is experimental: IAE(discovered-controller)/IAE(oracle) =
1.000 with settling 1.946 vs 1.945 s (exp11).

## P11 — Time-stencil sign discipline (engineering note)

A central difference stencil must be verified against a known function before use:
exp10 initially used a sign-flipped 4th-order time stencil, producing an exactly
negated Burgers equation (anti-dissipative; the rollout solver then diverged — itself
a useful built-in consistency check). Both stencils in this repo are now unit-tested
against sin(ωt). Lesson generalized: *in discovery pipelines, derivative operators are
part of the hypothesis class and deserve tests.*

## P12 — Why the integral form beats differentiation under state noise (sketch)

Let states be observed with iid noise ξ, Var(ξ)=σ². Differentiation amplifies:
Var(dξ/dt) ≈ 2σ²/dt² for central differences — at dt=2e-3 that multiplies noise
variance by 5e5 before regression. The window-integral formulation regresses
X(t+W)−X(t) against ∫Θ dt; differencing two noisy samples keeps variance 2σ²
(independent of dt), while each integrated feature column averages W samples, whose
noise contribution to the *feature* side scales like σ²·dt²·W/3-ish and, more
importantly, enters as an error-in-variables with mean zero uncorrelated across rows.
Net effect measured on LV: recovery rate at σ=0.05 rises from 7% to 63% and median
rollout error falls 12× (exp13). The trade-off: overlapping-window rows share cumulative
sums, so residuals are strongly correlated and naive significance thresholds are
optimistic; we therefore report the preregistered rollout endpoint, not p-values,
for this experiment. Full error analysis is left as future work; what is claimed here
is the empirical rescue plus the variance-scaling argument above.



## Honest-failure criterion (definition, not theorem)

A discovery attempt on dataset D is labeled **"no compact law found"** iff the best
candidate's held-out NMSE exceeds 1e-2 (preregistered). Rationale: a model that loses
>1% of variance out-of-sample has failed the predeclared usefulness bar regardless of
training fit; reporting it as a law would be curve fitting, not discovery. Applied to
adversarial cases A1/A2 in exp08.

