# T05.1 — fechamento de métricas, cobertura, artefatos e documentação

## Objetivo, base e núcleo protegido

T05.1 fecha a validação de trímeros sem mudar ciência. A base é `386772321537028a45c0958980b739edd3c5a780` em `main`, com 89 testes. Solver, equações, forças A/B/C, oráculo, regressões, geometrias e arquivos protegidos permanecem inalterados.

## Métricas

\[
F_{\mathrm{scale}}=\max_i\left(|\mathbf F_i^{\mathrm{ref}}|,|\mathbf F_i^{\mathrm{mod}}|ight),\qquad
F_{\mathrm{tol}}=128\,\epsilon_{\mathrm{mach}}F_{\mathrm{scale}},\quad
\epsilon_{\mathrm{mach}}=\operatorname{eps}(	exttt{float}).
\]

Não há piso absoluto `1.0`; vetor nulo satisfaz \(|\mathbf F|\le F_{\mathrm{tol}}\). O erro por partícula é

\[
arepsilon_{i,\mathrm{sym}}=rac{2|\mathbf F_i^{\mathrm{ref}}-\mathbf F_i^{\mathrm{mod}}|}{|\mathbf F_i^{\mathrm{ref}}|+|\mathbf F_i^{\mathrm{mod}}|}.
\]

Dois vetores nulos retornam zero; se qualquer vetor angular for nulo, o ângulo é `NaN`. A amplitude das correções é \(F_{\mathrm{RMS}}=[N^{-1}\sum_i|\mathbf F_i|^2]^{1/2}\), distinta de \(arepsilon_{\mathrm{RMS}}=[\sum_i|\Delta\mathbf F_i|^2/\sum_i|\mathbf F_i^{\mathrm{ref}}|^2]^{1/2}\).

## Validador, testes e simetrias

O validador usa RMS vetorial nas duas correções. A cobertura inclui vetores idênticos, opostos, ortogonais, nulos, resíduos globais, \(10^{-20}\), RMS, permutações, escalamento \(\lambda^2\), \(E_0\), \(f_0\), \(f_1\), rejeição de \(L_{\max}
e1\), cadeia e equilátero. A cadeia central tem erro zero e ângulo `NaN`; o equilátero preserva radialidade, soma nula e rotação de \(120^\circ\). Não se impõe soma nula ao escaleno.

## Artefatos, arquivos, verificações e limites

O sweep de 1.920 configurações foi regenerado; as duas RMS mudam por \(\sqrt2\), sem mudar A/B/C. O CSV de regressão tem nove linhas e hash `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`; T03/T04 preservam `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a` e `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`. Arquivos permitidos são métricas, exportação, teste, validador, artefatos T05 e documentação; arquivos científicos protegidos são proibidos. Aceite: testes sem warnings, hashes, `git diff --check`, identidade telescópica e sem Modelo D, \(N>3\) ou \(L_{\max}>1\). O limite científico permanece Rayleigh, \(N=3\), A/B/C e \(L_{\max}=1\).
