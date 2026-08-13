# Conventions

The physical pressure is represented by

\[
p(\mathbf r,t)=\operatorname{Re}\{p(\mathbf r)e^{-i\omega t}\}.
\]

For complex velocity amplitude \(v_0\), the reference energy density is
\(E_0=\rho_0|v_0|^2/4\).  It is always an explicit SI argument
(`J m^-3`) in the force API; no normalization is hidden in a global constant.

The positions are Cartesian 2D coordinates (metres) in a pressure **nodal**
plane.  For the force on probe `i` due to source `j`,
\(\mathbf d_{ij}=\mathbf r_i-\mathbf r_j\), so
\(\widehat{\mathbf d}_{ij}\) points from source to probe.  A negative signed
radial force is consequently attractive.  The implementation is valid only
for non-overlapping spheres, \(d\ge2a\).

The independent geometric parametrization is \(ka\) plus positions scaled by
\(a\); therefore \(kd=(ka)(d/a)\).  These three quantities must not be set as
independent parameters in sweeps.

This implementation follows the nodal pair interaction of Silva and Bruus,
*Physical Review E* 90, 063007 (2014), with the corrected sign audit and
two-body benchmark context recorded in Silva, *Brazilian Journal of Physics*
(2026), DOI: 10.1007/s13538-026-02102-x.  It does not apply to an antinodal
plane: that case has different incident-field physics and is not represented
by this API.

T02 also implements the corrected fifth-order analytical two-particle formula of the 2026 reference. Its signed radial component follows the same source-to-probe direction and non-overlap condition. Figure 2 uses \(100|F^{\mathrm{corr}}-F^{\mathrm{SB}}|/|F^{\mathrm{corr}}|\); it is undefined where the corrected force is zero.

## T03: Rayleigh multipolar solver

T03 provides a dense coupled solver at `Lmax=1`, with four modes and therefore `4N` complex field coefficients for `N` particles. The operator is oriented target <- source: rows are target modes, columns are source modes, and `R = source_position - target_position`. It solves `A = I - D_g U` with `numpy.linalg.solve`, thereby resumming all permitted rescattering orders. It produces field coefficients and diagnostics only; no multibody force is implemented.

SciPy provides spherical functions and complex Condon--Shortley harmonics; SymPy provides cached 3j-based Gaunt coefficients. The solver accepts nodal-plane centers only, while the low-level translation API is fully three-dimensional for reexpansion validation.

## T04 nodal interaction force

T04 implements the Model C Rayleigh interaction force specialized to the nodal plane: the external--scattered cross terms of Eqs. (22)/(27), not the unrestricted total force of Eq. (21), and without scattered--scattered quadratic terms. The T03 solver remains `Lmax_scatter=1`; reexpansion to `Lmax_evaluation=2` supplies only local regular field coefficients for force evaluation and never feeds back into the solver.

For each target, the self field is excluded. The documented Cartesian combinations preserve the derivation conjugation:

\[
F_x=\frac{\sqrt{30\pi}}{15}\,k a^3E_0
\operatorname{Re}\!\left[
f_1^*(b_{2,-1}-b_{2,1})
\right],
\]

\[
F_y=\frac{\sqrt{30\pi}}{15}\,k a^3E_0
\operatorname{Re}\!\left[
-i f_1^*(b_{2,1}+b_{2,-1})
\right].
\]

The API currently accepts real scalar \(f_1\), while production keeps `np.conj(f1)`.

## T05 trimer Model A/B/C comparison

For a cluster, Model A is the sum of T01 Silva--Bruus pair forces. Model B is the sum of isolated two-particle T04 solves, one solve for each unordered pair; it explicitly does not use the T02 truncated analytical formula. Model C is one global T04 solve of all three particles. The signed vector corrections are `Delta F^(2)=B-A` and `Delta F^(3)=C-B`, so `C-A=Delta F^(2)+Delta F^(3)` component by component.


Trimer reporting uses normalized `F/(a^2 E0)`. For one complete configuration,

\[
F_{\mathrm{scale}}=\max_i\left(|\mathbf F_i^{\mathrm{ref}}|,|\mathbf F_i^{\mathrm{mod}}|\right),\qquad
F_{\mathrm{tol}}=128\,\epsilon_{\mathrm{mach}}F_{\mathrm{scale}}.
\]

There is no absolute floor. A vector is numerically null when \( |\mathbf F|\le F_{\mathrm{tol}} \). Symmetric error is zero when both vectors are numerically null; angular error is `NaN` when either vector is null. Correction amplitude is \(F_{\mathrm{RMS}}=\sqrt{N^{-1}\sum_i|\mathbf F_i|^2}\), distinct from the dimensionless relative error \(\varepsilon_{\mathrm{RMS}}\). Model C minus B isolates collective rescattering inside the common Rayleigh \(L_{\max}=1\) basis; T05.1 changed derived metrics and artifacts, not A, B, or C.

## T06 connected quartet expansion

For a subset \(S\) containing particle \(i\), \(\mathbf F_i^C(S)\) denotes the Model C force obtained by solving only that subset with the same Rayleigh basis. Pair terms are \(\boldsymbol{\Phi}_i^{(2)}(\{i,j\})=\mathbf F_i^C(\{i,j\})\). For each embedded triplet,

\[
\boldsymbol{\Phi}_i^{(3)}(T)=\mathbf F_i^C(T)-\sum_{j\in T,\,j\ne i}\mathbf F_i^C(\{i,j\}).
\]

The reconstruction through three-body order and the recursive four-body term are

\[
\mathbf F_i^{(\le3)}=\mathbf F_i^B+\boldsymbol{\Phi}_{i,\Sigma}^{(3)},
\qquad
\boldsymbol{\Phi}_i^{(4)}=\mathbf F_i^C(Q)-\mathbf F_i^{(\le3)}.
\]

The equivalent closed form is

\[
\boldsymbol{\Phi}_i^{(4)}=\mathbf F_i^C(Q)
-\sum_{T\ni i,\,|T|=3}\mathbf F_i^C(T)
+\sum_{j\ne i}\mathbf F_i^C(\{i,j\}).
\]

This is a signed vector decomposition, not a decomposition of magnitudes. Dimensional amplitudes use \(F_{\mathrm{RMS}}=\sqrt{N^{-1}\sum_i|\mathbf F_i|^2}\), while \(\varepsilon_{\mathrm{RMS}}\) is a relative error. Numerical nullity and undefined angles continue to use the global T05.1 threshold without an absolute floor.

## T06.1 coupling and log-space scaling diagnostics

The minimum-distance coupling predictor and the collective geometric predictor
are

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3,
\qquad
\Lambda_i=|f_1|\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3,
\qquad
\Lambda_{\max}=\max_i\Lambda_i.
\]

Both are dimensionless. The geometric sum uses all neighbors and is evaluated
from Euclidean center distances. Non-overlap requires \(r_{ij}\ge 2a\).
For each fixed-shape, uniformly dilated T06 family,
\(\Lambda_{\max}=C_g\eta\); therefore replacing \(\eta\) by
\(\Lambda_{\max}\) changes only the intercept of an intrageometry power-law
fit, not its exponent or log-space residual diagnostics.

T06.1 fits strictly positive data without weighting according to

\[
\ln y=\ln C+p\ln x.
\]

For observed values \(y_n\) and fitted log values
\(\widehat{\ln y_n}\), the reported diagnostics are

\[
R^2_{\log}=1-
\frac{\sum_n(\ln y_n-\widehat{\ln y_n})^2}
{\sum_n(\ln y_n-\overline{\ln y})^2},
\]

\[
\operatorname{RMSE}_{\log}=
\sqrt{\frac{1}{N}\sum_n(\ln y_n-\widehat{\ln y_n})^2}.
\]

The exponent is \(p\) when the predictor is \(x\). When the horizontal axis is
\(x^2\), the displayed exponent relative to that transformed predictor is
\(q=p/2\). These diagnostics describe amplitude-ratio collapse; they do not
establish an additive force fraction or a universal validity threshold.

## T07 multipolar Model D

For \(\ell\geq1\), Model D uses the leading Rayleigh coefficient

\[
s_\ell=i\frac{3\ell f_1}{(2\ell-1)!!(2\ell+1)!![2(2\ell+1)-(\ell-1)f_1]}(ka)^{2\ell+1},
\]

while \(s_0=-if_0(ka)^3/3\). The complete ordering contains \((L_{\max}+1)^2\) modes per particle. Planar reflection permits \(\ell+m\) odd; inactive modes are returned as exact zeros. This is a symmetry of the nodal field, not a rule that only odd \(\ell\) may occur in a generic planar cluster.

The physical and balanced systems are

\[
(\mathbf I-\mathbf D\mathbf U)\mathbf s=\mathbf D\mathbf a_{\rm ext},
\qquad
(\mathbf I-\mathbf D^{1/2}\mathbf U\mathbf D^{1/2})\mathbf q=\mathbf D^{1/2}\mathbf a_{\rm ext}.
\]

The physical residual is measured in the first equation, while reported conditioning refers primarily to the balanced matrix. Multipole order does not count re-scattering events: the coupled solve already resums all paths admitted at fixed \(L_{\max}\).

Connected terms at order \(L\) are defined only from Model-D subset solutions at that same \(L\), by vector inclusion--exclusion. Successive convergence is measured separately for the total force and each nonzero connected term. A relative ratio is undefined for a numerically null term; its applicability flag must be false rather than silently interpreting zero as a measured ratio.

## T08 transferability conventions

The historical Model B remains the sum of isolated Rayleigh dimers at
\(L=1\). T08 introduces a separate, explicitly order-matched diagnostic,

\[
\mathbf F_i^{B_L}=\sum_{j\ne i}\mathbf F_{ij}^{D,N=2,L},
\]

called the multipolarly matched pairwise baseline. Thus \(B_1=B\), while
\(B_L=D_L\) for a dimer. It does not redefine Model B. The signed vector
identity audited in every configuration is

\[
\mathbf F^{D_L}-\mathbf F^A=
(\mathbf F^{B_L}-\mathbf F^A)+(\mathbf F^{D_L}-\mathbf F^{B_L}).
\]

With the project RMS vector magnitude, the primary and residual errors are

\[
\varepsilon_A=\frac{F_{\mathrm{RMS}}(\mathbf F^A-\mathbf F^D)}
{F_{\mathrm{RMS}}(\mathbf F^D)},
\qquad
\varepsilon_B=\frac{F_{\mathrm{RMS}}(\mathbf F^{B_L}-\mathbf F^D)}
{F_{\mathrm{RMS}}(\mathbf F^D)}.
\]

The associated amplitudes are

\[
Y_{\mathrm{2B}}=\frac{F_{\mathrm{RMS}}(\mathbf F^{B_L}-\mathbf F^A)}
{F_{\mathrm{RMS}}(\mathbf F^D)},\quad
Y_{\mathrm{coll}}=\frac{F_{\mathrm{RMS}}(\mathbf F^D-\mathbf F^{B_L})}
{F_{\mathrm{RMS}}(\mathbf F^D)},\quad
Y_{\mathrm{mp}}=\frac{F_{\mathrm{RMS}}(\mathbf F^D-\mathbf F^{D_1})}
{F_{\mathrm{RMS}}(\mathbf F^D)}.
\]

The three dimensionless predictors are

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3,
\qquad
\Lambda_{\max}=|f_1|\max_i\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3,
\]

\[
\rho_1=\rho(\mathbf K_b^{(1)}),
\qquad
\mathbf K_b^{(1)}=\mathbf I-\mathbf A_b^{(1)}.
\]

Here \(\rho_1\) is computed from the balanced \(L=1\) rescattering operator,
not from the system matrix itself or from its condition number.

Convergence is assessed independently for \(D_L\), \(B_L\), and
\(R_L=D_L-B_L\). A quantity is confirmed only when its last two successive
normalized RMS changes are both at most \(10^{-3}\). The collective residual
is numerically resolved only when \(\varepsilon_B>5u_R\), with \(u_R\) the
larger of its last two truncation changes. Non-applicable ratios use an
explicit flag and a finite placeholder; they are never represented by
`NaN` or `inf`.

Calibration is restricted to \(N\leq4\). Clusters with \(N=6,10\) form an
external holdout and cannot select a predictor, fit a power law, or define a
threshold. The thresholds at 1%, 5%, and 10% are conservative empirical
nodal-plane thresholds within the sampled domain, not universal constants.

## T09 analytical rho_1 conventions

For \(f_0=0\), \(L=1\), and centers in the nodal plane, the only active
scattered channel per particle is \((\ell,m)=(1,0)\). The balanced operator
therefore reduces exactly, within the leading-Rayleigh dipole model, to

\[
(K_b)_{ii}=0,\qquad
(K_b)_{ij}=
\frac{f_1}{2}\left(\frac{a}{r_{ij}}\right)^3
e^{ikr_{ij}}(1-ikr_{ij}).
\]

The word `exact` in this section means algebraically identical to the
implemented \(L=1\) leading-Rayleigh operator. It does not mean an exact Mie
coefficient or complete radiation force.

The near-field operator drops only the retarded factor:

\[
(K_b^{\mathrm{nf}})_{ij}
=
\frac{f_1}{2}\left(\frac{a}{r_{ij}}\right)^3.
\]

The production value of \(\rho_1\) always uses the retarded operator, not this
near-field approximation. The Neumann index \(p\) counts additional
rescattering events; it is unrelated to multipole truncation \(L_{\max}\).
Spectral-radius convergence and force accuracy must also remain distinct:
\(\rho_1<1\) guarantees convergence of the finite-dimensional matrix Neumann
series for every source, not a prescribed Silva--Bruus error.

## T10 isolated-sphere Mie conventions

The exact lossless-fluid coefficient uses the same \(e^{-i\omega t}\)
convention and outgoing Hankel function as the rest of the project. Define

\[
x=ka,\quad y=x\sqrt{\widetilde\rho\widetilde\kappa},\quad
\beta=\sqrt{\widetilde\kappa/\widetilde\rho}.
\]

Then

\[
s_\ell=-\frac{\beta j_\ell(x)j_\ell'(y)-j_\ell(y)j_\ell'(x)}
{\beta h_\ell^{(1)}(x)j_\ell'(y)-j_\ell(y){h_\ell^{(1)}}'(x)}.
\]

Material ratios obey
\(\widetilde\kappa=1-f_0\),
\(\widetilde\rho=(2+f_1)/[2(1-f_1)]\), and
\(c_p/c_0=(\widetilde\rho\widetilde\kappa)^{-1/2}\) for
\(-2<f_1<1\). Exactly \(f_1=1\), without clipping, selects
\(s_\ell=-j_\ell'(x)/{h_\ell^{(1)}}'(x)\). Undefined tabulated quantities
use `NaN` plus an explicit applicability flag. The exact-Mie label applies
only to the isolated-sphere T-matrix and not to the collective force.
## T11 complete-reference Model E

For each particle, \(a\) denotes the external incident BSCs, \(b\) the
effective incident BSCs, and \(d=Db\) the scattered BSCs. The legacy
effective-incident system is

\[
A_b b=a,
\qquad A_b=I-UD.
\]

The equivalent scattered-field system is

\[
A_d d=Da,
\qquad A_d=I-DU.
\]

Using the elementwise principal complex square root,

\[
S=D^{1/2},
\qquad
A_q=I-SUS,
\qquad
A_q q=Sa.
\]

Production solves only the balanced system with `numpy.linalg.solve` and
reconstructs without division by \(S\):

\[
d=Sq,
\qquad
b=a+Ud.
\]

No inverse, pseudoinverse, least-squares solve, or magnitude pruning is used.
Legacy public attributes continue to denote \(A_b\), \(a\),
\(\kappa(A_b)\), and the residual of the legacy equation.

The balanced backward error is

\[
\eta_q=
\frac{\lVert A_q q-Sa\rVert}
{\lVert A_q\rVert\lVert q\rVert+\lVert Sa\rVert}.
\]

The physical closure diagnostics are

\[
r_b=\frac{\lVert b-a-Ud\rVert}
{\lVert b\rVert+\lVert a\rVert+\lVert Ud\rVert},
\qquad

For \(n=0,\ldots,L_{\max}-1\), define

\[
\Gamma_n=s_n+s_{n+1}^*+2s_ns_{n+1}^*.
\]

The complete force is

\[
F_x+iF_y=\frac{iE_0}{k^2}\sum_{n,m}
\sqrt{\frac{(n+m+1)(n+m+2)}{(2n+1)(2n+3)}}
\left[\Gamma_n b_{nm}b_{n+1,m+1}^*
+\Gamma_n^*b_{n,-m}^*b_{n+1,-m-1}\right],
\]

\[
F_z=\frac{2E_0}{k^2}\operatorname{Im}\sum_{n,m}
\sqrt{\frac{(n-m+1)(n+m+1)}{(2n+1)(2n+3)}}
\Gamma_n b_{nm}b_{n+1,m}^*.
\]

These prefactors already use \(E_{\mathrm{LAS}}=2E_0\). With \(c=b-a\),

\[
F_{\mathrm{total}}=\mathcal F[b],\quad
F_{\mathrm{external}}=\mathcal F[a],\quad
F_{\mathrm{int}}=F_{\mathrm{total}}-F_{\mathrm{external}},
\]

\[
F_{\mathrm{ss}}=\mathcal F[c],\qquad
F_{\mathrm{ext-sc}}=F_{\mathrm{int}}-F_{\mathrm{ss}}.
\]

The quadratic recoil \(2s_ns_{n+1}^*\) inside \(\Gamma_n\) is not the same
object as \(\mathcal F[c]\). The principal Model-E API requires
\(L_{\max}\ge2\), and every force uses coefficients through precisely the
requested truncation; coefficients of order \(L_{\max}+1\) are never
invented. Convergence requires two consecutive applicable changes not
exceeding \(10^{-5}\), separately for each force channel.

## P1.2 complete-dimer Model-E baseline

The new editorial baseline is

\[
\mathbf F_i^{B_E}
=
\sum_{j\ne i}\mathbf F_i^{E,\{i,j\}},
\]

where each unordered pair \((i,j)\), \(i<j\), is solved as an isolated
two-particle Model-E problem. Historical Model B remains the isolated
Rayleigh \(L=1\) sum and \(B_L\) remains the order-matched leading-Rayleigh
diagnostic; neither name nor implementation is changed.

`solve_model_be_nodal` preserves the original particle indices and passes
the two original Cartesian positions to Model E in the order \((i,j)\). Pair
records are emitted in lexicographic \((i,j)\) order. The force contributed to
\(B_E\) is exactly the pair's final Model-E `interaction_forces_xyz`, with
its first row accumulated on particle \(i\) and its second row on particle
\(j\). No coordinate rotation, sign reconstruction or action--reaction
assumption is used by the aggregator.

Each pair evaluates its own integer sequence
\(L_{\max}=2,\ldots,21\). It cannot stop before \(L_{\max}=5\) and stops only
after total, interaction, external--scattered and scattered--scattered force
channels have each closed two consecutive applicable normalized RMS changes
at or below \(10^{-5}\). A numerically null channel retains
`applicable=false` and does not fabricate a relative error.

The final pair result must also pass the established Model-E gates:
`balanced_sqrt`, finite force/coefficient/diagnostic values, balanced
condition number below 10, balanced backward error and the three
closure/decomposition residuals below \(10^{-12}\), consistent mode
dimension, and the scale-relative planar \(F_z\) tolerance. These gates are
implemented once in `evaluate_model_e_numerical_diagnostics` and reused by
\(B_E\).

`ModelBEResult.pair_ledger` retains every attempted pair, its two individual
dimer forces, attempted/evaluated/final/failed orders, channel histories,
applicability, convergence, diagnostics and explicit failure stage/reason.
A failure is local to its pair and later pairs are still audited. Global
`eligible` is true only if every pair is eligible; otherwise
`forces_xyz=None`, so no partial or imputed \(B_E\) vector can escape.

P1.2 defines and unit-tests this orchestration contract with an injected fake
solver. Common-order sensitivity, full physical dimer identities,
rotation/reflection, action--reaction and asymptotic-limit evidence remain
assigned to P1.3.

## T12 Model-E sentinel comparisons

The T12 reference is the interaction force

\[
\mathbf F^E=\mathbf F^E_{\mathrm{int}}
=\mathbf F^E_{\mathrm{ext-sc}}+\mathbf F^E_{\mathrm{ss}},
\]

not the total force. Planar A and D vectors are padded with an explicit zero
\(z\) component and are compared with complete three-dimensional E vectors.
For any vector field,

\[
\mathcal R(\mathbf F)=
\left[\frac1N\sum_i\lVert\mathbf F_i\rVert_2^2\right]^{1/2}.
\]

The principal error is

\[
\varepsilon_A^E=
\frac{\mathcal R(\mathbf F^A-\mathbf F^E)}{\mathcal R(\mathbf F^E)}.
\]

The exact signed-vector mechanism identity is

\[
\mathbf F^E-\mathbf F^A=
(\mathbf F^D-\mathbf F^A)
+(\mathbf F^E_{\mathrm{ext-sc}}-\mathbf F^D)
+\mathbf F^E_{\mathrm{ss}}.
\]

Normalized RMS amplitudes of these terms are diagnostics, not additive
fractions. A ratio is applicable only when its reference scale exceeds the
global scale-relative tolerance \(128\epsilon_{\mathrm{mach}}F_{\mathrm{scale}}\).
Otherwise its stored value is zero with an explicit flag and reason; no
absolute floor is used.

## T12.1 signed-mechanism and predictor diagnostics

Define the signed vector fields

\[
\mathbf C_D=\mathbf F^D-\mathbf F^A,\qquad
\mathbf C_M=\mathbf F^E_{\mathrm{ext-sc}}-\mathbf F^D,\qquad
\mathbf C_S=\mathbf F^E_{\mathrm{ss}},
\]

\[
\mathbf C=\mathbf F^E-\mathbf F^A
=\mathbf C_D+\mathbf C_M+\mathbf C_S.
\]

Their vector-field inner product is

\[
\langle\mathbf X,\mathbf Y\rangle
=\frac1N\sum_i\mathbf X_i\mathbin{\cdot}\mathbf Y_i.
\]

Cosines use
\(\mu_{XY}=\langle X,Y\rangle/[\mathcal R(X)\mathcal R(Y)]\), and signed
projections use
\(p_X=\langle X,C\rangle/\langle C,C\rangle\), so
\(p_D+p_M+p_S=1\) when applicable. Amplitude ratios such as
\(\mathcal R(C_S)/\mathcal R(C_D)\) use the same scale-relative numerical
nullity convention as T12 and never an absolute floor.

Every fitted candidate uses unweighted least squares in logarithmic space,

\[
\ln y=\ln C+p\ln x.
\]

Predictive diagnostics are out-of-fold under deterministic
leave-\((N,\mathrm{family})\)-out validation. P0 is the unchanged frozen T08
law; P1, P2, P3, and P4 use respectively \(\eta\),
\(\Lambda_{\max}\), \(\rho_1\), and \(\varepsilon_A^D\). P4 is explicitly a
reference-derived diagnostic rather than a standalone validity predictor.

## T13 external-validation conventions

The external predictor and its conservative envelope are frozen as

\[
\widehat\varepsilon_{M1}
=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\qquad
\widehat\varepsilon_{M1,\mathrm{safe}}
=2.5699703122019222\,\widehat\varepsilon_{M1}.
\]

P3 is reported only as a frozen comparator and never participates in the M1
gate. Predicted and observed safety both use a strict inequality against the
tolerance; equality is unsafe. The external response is exclusively
\(\varepsilon_A^E\) formed from the complete Model-E interaction force.

An external case is eligible exactly when the interaction channel has two
successive applicable changes no larger than \(10^{-5}\), every numerical
diagnostic passes, and the interaction-force error is applicable. A case that
fails any item remains in the 24-row summary but is excluded from metrics and
threshold audits without imputation. The standard order cap is 13 and the
interaction-only extension cap is 21.

Multiplicative metrics use positive errors and predictions in natural-log
space. Conservative false-safe means
\(\widehat\varepsilon_{\mathrm{safe}}<\tau\) while
\(\varepsilon_A^E\geq\tau\). Mechanism amplitudes remain RMS vector-field
diagnostics and are not additive scalar force fractions.

## T14 scale-out conventions

The scale-out geometry uses the direct geometric coupling

\[
\Lambda_{\max}=|f_1|a^3\max_i\sum_{j\ne i}r_{ij}^{-3}.
\]

For each centered unit-minimum template, define
\(S_{N,g}=\max_i\sum_{j\ne i}r_{ij}^{-3}\) and apply the analytic scale

\[
d=\left(\frac{|f_1|S_{N,g}}
{\Lambda_{\max}^{\mathrm{target}}}\right)^{1/3}.
\]

This changes particle count and geometric aperture while preserving the
specified coupling target. The confirmatory response remains
\(\varepsilon_A^E=\mathcal R(\mathbf F^A-\mathbf F^E_{\mathrm{int}})/
\mathcal R(\mathbf F^E_{\mathrm{int}})\), in three dimensions. A case enters
metrics only when the interaction channel has two-step convergence, all
numerical diagnostics pass, and the error denominator is applicable.

The frozen M1 and diagnostic P3 laws, safety factors, strict inequalities,
and 1%, 5%, and 10% tolerances are identical to T13. P3 cannot affect the M1
gate. T14 uses orders \(L=2,\ldots,13\), stops from \(L=5\) only after all
four force channels are confirmed, and never imputes an unconfirmed force.

## T14.1 large-\(N\) confirmation conventions

T14.1 retains the frozen T12.3--T14 response and predictor definitions. For
each particle,

\[
\Lambda_i=|f_1|\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3,
\qquad
\Lambda_{\max}=\max_i\Lambda_i,
\]

and the confirmatory response is

\[
\varepsilon_A^E=
\frac{\operatorname{RMS}(\mathbf F^A-\mathbf F^E_{\mathrm{int}})}
{\operatorname{RMS}(\mathbf F^E_{\mathrm{int}})}.
\]

The complete ordered vector \(\{\Lambda_i\}\), its SHA-256, minimum, mean,
median, standard deviation, percentiles 10 and 90, maximum,
\(\overline\Lambda/\Lambda_{\max}\), fraction satisfying
\(\Lambda_i\ge0.9\Lambda_{\max}\), and first maximizing index are frozen
before Model-E evaluation. Each fixed template is uniformly scaled to one of
the same four T13--T14 targets.

Convergence uses every integer order \(L_{\max}=2,3,\ldots,13\), never stops
before order 5, and requires two successive applicable changes no larger than
\(10^{-5}\) for total, interaction, external--scattered, and
scattered--scattered channels. The order-level numerical gate uses the
`balanced_sqrt` solver, balanced condition number below 10, balanced backward
error and all three closure/decomposition errors below \(10^{-12}\), finite
outputs, consistent mode dimension, and planar symmetry. The separately
recorded physical residual of the unbalanced equation is diagnostic and is
not substituted for this preregistered balanced gate.

A case is eligible only after the campaign and interaction convergence are
confirmed, the final numerical diagnostics pass, the A--E denominator is
applicable, and every phase-A coordinate, coupling, prediction, protocol, and
hash identity remains intact. No missing response is imputed. M1 alone is
confirmatory; P3 remains diagnostic and cannot change the gate.
