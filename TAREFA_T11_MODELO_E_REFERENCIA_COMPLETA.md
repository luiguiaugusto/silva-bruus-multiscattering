# T11 — Modelo E e referência multipolar completa

## Objetivo e base

A T11 parte do commit \`153403c59571ad081098248165b5b184c3721179\`
e implementa uma referência coletiva completa para esferas fluidas idênticas,
sem alterar os Modelos A--D. A linha de base foi de 258 testes.

\[
\boxed{\text{Modelo E}=\text{Mie exato}+\text{múltiplo espalhamento global}+\text{força multipolar completa}.}
\]

O escopo permanece restrito a fluido ideal, esferas sem perdas, fixas, não
sobrepostas e centradas no plano nodal.

## Sistema global e T-matrix

Com \(a\) externo, \(b\) incidente efetivo, \(d\) espalhado e a T-matrix
diagonal exata \(D\) da T10,

\[
b=a+Ud,\qquad d=Db,\qquad \boxed{(I-UD)b=a}.
\]

A produção usa \`numpy.linalg.solve\` e nunca forma uma inversa. Um teste
independente resolve \((I-DU)d=Da\). A tradução mantém a orientação
\`target <- source\`, e a T-matrix vem exclusivamente de
\`mie_scattering_coefficients_from_contrasts\`.

A ordenação completa é \(n^2+n+m\). A redução planar conserva \(n+m\) ímpar,
retorna vetores na base completa e preenche canais inativos com zero. Ela é
comparada com a solução de base completa. Nenhum modo é podado por limiar
numérico.

## Funcional de força

Para \(n=0,\ldots,M-1\),

\[
\Gamma_n=s_n+s_{n+1}^*+2s_ns_{n+1}^*,
\]

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

Os prefatores incorporam \(E_{\mathrm{LAS}}=2E_0\). A API principal exige
\(L_{\max}\ge2\), pois a força acopla ordens adjacentes. Nenhum coeficiente de
ordem \(L_{\max}+1\) é inventado.

## Decomposição dos canais

Com \(c=b-a\),

\[
F_{\mathrm{total}}=\mathcal F[b],\qquad
F_{\mathrm{external}}=\mathcal F[a],\qquad
F_{\mathrm{int}}=F_{\mathrm{total}}-F_{\mathrm{external}},
\]

\[
F_{\mathrm{ss}}=\mathcal F[c],\qquad
F_{\mathrm{ext-sc}}=F_{\mathrm{int}}-F_{\mathrm{ss}}.
\]

Portanto \(F_{\mathrm{int}}=F_{\mathrm{ext-sc}}+F_{\mathrm{ss}}\). O termo
\(\mathcal F[c]\) é distinto do recuo \(2s_ns_{n+1}^*\) dentro de
\(\Gamma_n\).

## Oráculo independente e redução Rayleigh

O oráculo de validação reconstrói \(\psi\) e \(g=k^{-1}\nabla\psi\) e integra

\[
\frac{\overline S}{E_0}=-(|g|^2-|\psi|^2)I+2\operatorname{Re}(gg^\dagger),
\qquad
F=-E_0R^2\int(\overline S/E_0)\cdot\widehat e_r\,d\Omega.
\]

Foram usados \(R/a=1.01,1.04\) e quadraturas \(24\times48\) e
\(32\times64\). O maior erro relativo de componente resolvida foi
\(6.30\times10^{-15}\). O teste artificial \(s_0=s_2=0\), com \(s_1\) de
Rayleigh, reproduziu a força external--scattered do Modelo D em \(L=1\).

## Campanha compacta

A campanha contém seis casos: dímero axial com \(d/a=2.5\), dímero diagonal
com \(d/a=4\), dímero rígido com \(d/a=3\), equilátero de lado \(3a\),
escaleno com \(d_{\min}/a=2.7\) e quarteto irregular com
\(d_{\min}/a=2.8\). O dímero diagonal usa \(ka=0.05,f_1=0.4\); o rígido usa
\(ka=0.1,f_1=1\); os demais usam \(ka=0.1,f_1=0.8\). Todos usam
\(f_0=0\) e \(a=E_0=1\).

As ordens \(L_{\max}=2,\ldots,7\) foram obrigatórias. Todos os casos foram
estendidos até \(L_{\max}=9\) porque pelo menos um canal ainda não tinha duas
confirmações. A convergência exige duas mudanças sucessivas aplicáveis
\(u_L\le10^{-5}\), separadamente para total, interação, external--scattered e
scattered--scattered.

O dímero diagonal confirmou os quatro canais em \(L=7,7,6,9\). Dímero rígido,
equilátero e quarteto irregular confirmaram total/interação em \(L=9\) e
external--scattered em \(L=8\), mas não scattered--scattered. O escaleno
confirmou apenas external--scattered, em \(L=9\). O dímero axial não confirmou
nenhum canal até o limite. Casos assim são rotulados não confirmados, nunca
divergentes.

O maior resíduo foi \(9.58\times10^{-5}\), no dímero axial em \(L=9\), e o
maior número de condição foi \(5.91\times10^{24}\). Esses diagnósticos são
preservados sem clipping, regularização ou troca silenciosa do sistema.

## Testes e artefatos

Os testes cobrem identidades lineares, base completa/reduzida, transparência,
rigidez, partícula isolada, validação de entradas, simetrias, energia, fase,
rotação, translação, permutação, redução ao Modelo D, tensor de tensões,
decomposição, finitude e determinismo. A suíte final tem 295 testes sem
warnings.

\`\`\`text
f017993c893a6a1d8db5161007ba8361dc55770922849101e8ac193d45ccf893  results/data/t11_model_e_convergence.csv
d1a8e89c62a248ac339d5f8c1e51c35b30651034460ca5b2bed0419d05f585fb  results/data/t11_force_oracle.csv
6c2b83134b306ab6b113e058439348a9de21f5f735a468a891e0868ad6f53986  results/data/t11_force_decomposition.csv
b96f9e089a833b175bea0621fef6b708c66f73a47f523655300c24a67ac9301f  results/figures/t11_model_e_validation.png
\`\`\`

Duas execuções no mesmo ambiente foram byte-idênticas. A figura foi
inspecionada quanto a legibilidade, escalas, linha de tolerância, curvas e
ausência de \`NaN\` ou \`inf\` físico.

## Limitações e próximo passo

Não foram adicionados absorção, viscosidade, streaming, paredes, elasticidade,
esferas não idênticas ou movimento. A T11 não executa sentinelas T12, não
recalibra \(\rho_1\) e não abre o holdout T13--T14.

\[
\boxed{\text{convergência interna do Modelo E}\ne
\text{validação dos limiares de }\rho_1.}
\]
