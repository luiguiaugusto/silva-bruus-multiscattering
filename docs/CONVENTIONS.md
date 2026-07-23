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
