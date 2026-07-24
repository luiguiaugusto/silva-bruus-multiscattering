# T05.1 — fechamento de métricas, testes, artefatos e documentação

## Objetivo, escopo e núcleo protegido

A T05.1 fechou a validação da comparação de trímeros aprovada na T05, cujo commit-base é `386772321537028a45c0958980b739edd3c5a780`. Foram permitidos somente métricas derivadas, sua exportação, o teste multibody, o validador T05, o sweep e a figura T05, e documentação. O solver, equações, geometrias, forças A, B e C, oráculo escalar e regressões científicas ficaram protegidos.

## Métricas finais

Para a configuração completa,

\[
F_{\mathrm{scale}}=\max_i\left(
|\mathbf F_i^{\mathrm{ref}}|,
|\mathbf F_i^{\mathrm{mod}}|
\right),
\qquad
F_{\mathrm{tol}}=128\,\epsilon_{\mathrm{mach}}F_{\mathrm{scale}},
\]

com \(\epsilon_{\mathrm{mach}}=\operatorname{eps}(\texttt{float})\). Não há piso absoluto `1.0`; um vetor é numericamente nulo quando \(|\mathbf F|\le F_{\mathrm{tol}}\). O erro simétrico é

\[
\varepsilon_{i,\mathrm{sym}}=
\frac{2|\mathbf F_i^{\mathrm{ref}}-\mathbf F_i^{\mathrm{mod}}|}
{|\mathbf F_i^{\mathrm{ref}}|+|\mathbf F_i^{\mathrm{mod}}|}.
\]

Quando ambos os vetores são numericamente nulos ele é zero. Quando qualquer vetor necessário à comparação angular é nulo, o ângulo é `NaN`. Isto evita direções espúrias e não apaga o caso físico de forças opostas de ordem \(10^{-20}\). A amplitude RMS é

\[
F_{\mathrm{RMS}}=
\sqrt{\frac{1}{N}\sum_{i=1}^{N}|\mathbf F_i|^2},
\]

distinta do erro relativo

\[
\varepsilon_{\mathrm{RMS}}=
\left[
\frac{\sum_i|\mathbf F_i^{\mathrm{ref}}-\mathbf F_i^{\mathrm{mod}}|^2}
{\sum_i|\mathbf F_i^{\mathrm{ref}}|^2}
\right]^{1/2}.
\]

## Testes, simetrias e artefatos

Foram cobertos vetores idênticos, opostos, ortogonais e nulos; resíduos relativos à escala global; forças \(10^{-20}\); RMS vetorial; duas permutações; escalamento dimensional \(\lambda^2\); dependência em \(E_0\), \(f_0\) e \(f_1\); rejeição de \(L_{\max}\ne1\); e simetrias da cadeia e do equilátero. A cadeia central tem erro zero e ângulo `NaN`; o equilátero tem módulos iguais, radialidade, soma nula e rotação de \(120^\circ\). Não se impõe soma nula ao escaleno.

O validador passou a usar RMS vetorial nas duas correções. O sweep de 1.920 configurações foi auditado; as duas colunas RMS obedecem \(F_{\mathrm{RMS,new}}=\sqrt{2}F_{\mathrm{RMS,old}}\), enquanto A, B e C não mudam. Aceite: testes sem warnings, `git diff --check`, hashes protegidos, identidade telescópica e nenhuma extensão para Modelo D, \(N>3\) ou \(L_{\max}>1\). Os limites são regime de Rayleigh, \(N=3\) e somente Modelos A, B e C.
