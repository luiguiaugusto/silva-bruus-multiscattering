# T06.1 — análise de escala da expansão conectada

## Objetivo

A T06.1 pós-processa os 1.920 casos já versionados da T06 para avaliar a
hierarquia de escala das amplitudes conectadas de três e quatro corpos. Nenhuma
varredura de força é refeita. A única avaliação adicional do Modelo C é um
dímero centrado na separação \(d/a=2.1\), usado na tabela comparativa
\(N=2,3,4\).

## Escopo científico

O conjunto permanece em \(a=E_0=1\), \(ka=0.1\), \(f_0=0\),
\(f_1\in\{0.1,0.4,0.8,1.0\}\) e \(d_{\min}/a\in[2.1,10]\). São analisadas
somente as famílias `linear_chain`, `square` e `irregular`, todas planares e
uniformemente dilatadas, com Modelo C Rayleigh em
\(L_{\max}^{\mathrm{scatter}}=1\).

Os preditores são

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3,
\qquad
\Lambda_i=|f_1|\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3,
\qquad
\Lambda_{\max}=\max_i\Lambda_i.
\]

Para cada resposta positiva, o ajuste não ponderado usa todos os pontos:

\[
\ln y=\ln C+p\ln x.
\]

São registrados \(C\), \(p\), \(R^2_{\log}\), a RMSE logarítmica e o maior
resíduo logarítmico absoluto. Não há descarte, binning, suavização, ponderação,
ajuste robusto ou inferência estatística.

## Respostas e critérios

As respostas são

\[
Y_3=
\frac{F_{\mathrm{RMS}}(\boldsymbol{\Phi}_{\Sigma}^{(3)})}
{F_{\mathrm{RMS}}(\mathbf F^C)},
\qquad
Y_4=
\frac{F_{\mathrm{RMS}}(\boldsymbol{\Phi}^{(4)})}
{F_{\mathrm{RMS}}(\mathbf F^C)}.
\]

Devem ser produzidos 16 ajustes: quatro grupos, dois preditores e duas
respostas. Os fatores geométricos esperados
\(C_g=\Lambda_{\max}/\eta\) são \(2.125\), \(2.353553390593274\) e
\(1.996580257145743\) para cadeia, quadrado e quadrilátero irregular.

O aceite exige os expoentes auditados, as reduções agrupadas da RMSE
logarítmica, CSVs com 16, 2 e 7 linhas, duas figuras determinísticas, testes com
warnings tratados como erros e preservação dos nove hashes das T03–T06.

## Interpretação e limitações

Os dados sustentam descritivamente
\(Y_3\propto\eta^{0.94\text{–}0.95}\) e
\(Y_4\propto\eta^{1.91\text{–}1.92}\), consistentes com ordens aproximadas
\(O(\eta)\) e \(O(\eta^2)\). Como
\(\Lambda_{\max}=C_g\eta\) dentro de cada família fixa, a troca de preditor não
altera expoente ou qualidade do ajuste intrageometria; ela apenas pode reduzir
diferenças de intercepto no ajuste agrupado.

Essa melhora não define um critério universal. A análise cobre apenas
\(ka=0.1\), \(f_1>0\), três formas fixas, \(N\le4\), plano nodal, regime de
Rayleigh e força external–scattered. O acoplamento completo contém funções de
Hankel, \(kd\) varia, e não há demonstração analítica. A amplitude da soma
vetorial \(\boldsymbol{\Phi}_{\Sigma}^{(3)}\) não é a soma das amplitudes dos
trímeros. Modelo D, multipolos superiores, T07 e limiares de validade
permanecem fora do escopo.

## Implementação e artefatos

A API pública criada é composta por `coupling_eta`,
`maximum_geometric_coupling`, `fit_power_law` e a dataclass congelada
`PowerLawFit`. O script `scripts/analyze_t06_scaling.py` lê os quatro CSVs das
T05/T06, valida conteúdo, ordem e finitude, executa os 16 ajustes e escreve:

- `results/data/t06_1_scaling_fits.csv`, com 16 linhas;
- `results/data/t06_1_collapse_summary.csv`, com duas linhas;
- `results/data/t06_1_body_order_summary.csv`, com sete linhas;
- `results/figures/t06_1_eta_scaling.png`;
- `results/figures/t06_1_lambda_scaling.png`.

A suíte deve cobrir validação das entradas, invariâncias de
\(\Lambda_{\max}\), fatores geométricos canônicos, ajuste sintético exato,
integridade dos 1.920 casos e todas as regressões numéricas. Duas execuções no
mesmo ambiente devem gerar hashes idênticos.

## Integridade protegida

Os hashes oficiais que não podem mudar são:

- T03 CSV: `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`;
- T04 CSV: `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`;
- T05 regression: `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`;
- T05 sweep: `dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1`;
- T05 figure: `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`;
- T06 regression: `8d05db59dc4a44ee76118537af40db76aa386c8098f3f87f4359830cb5f9dea0`;
- T06 sweep: `36f64ebd16ea1df52bf4074d42bd83356306bfc17c613267b0def746901689b5`;
- T06 model figure: `547945fa0fd658565fb837416f8f4a5c65bd4963c0e2e50226510f82f7af17d0`;
- T06 body figure: `e4ce5f5a0c5f1d212d72ddd1bc6d35e1ecded8e86e437efeaf4a19bce4ab6d16`.

O aceite final exige testes sem warnings, `git diff --check`, ausência de
caracteres ASCII de controle nos documentos, inspeção visual das figuras,
escopo Git restrito aos arquivos permitidos, commit
`feat: analyze connected body-order scaling` e push para `main`.
