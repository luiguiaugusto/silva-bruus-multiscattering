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

The API currently accepts real scalar (f_1), while production keeps `np.conj(f1)`.
