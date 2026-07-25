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
