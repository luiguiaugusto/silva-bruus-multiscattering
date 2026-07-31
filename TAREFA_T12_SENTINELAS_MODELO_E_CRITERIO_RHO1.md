# T12 — sentinelas do Modelo E para o critério \(\rho_1\)

## Objetivo e proveniência

A T12 audita, sem recalibração, se a lei e os limiares empíricos de
\(\rho_1\) congelados pela T08 continuam diagnosticando de modo conservador o
erro de Silva–Bruus quando a referência passa do Modelo D para o Modelo E. A
base é o fechamento estabilizado da T11.1. O domínio da auditoria é somente a
calibração \(N\leq4\); nenhum caso \(N=6\) ou \(N=10\) é resolvido pelo Modelo
E nesta tarefa.

Os 28 casos foram pré-registrados antes de qualquer avaliação E: são quatro
faixas de \(\rho_1\) em cada um dos sete estratos
\((N,\text{família})\). O arquivo `t12_sentinel_manifest.csv` preserva a ordem,
os parâmetros, a ordem de referência D e a proveniência
`pre_registered_t08_calibration`.

## Núcleo científico protegido

A T12 não altera os Modelos A–E, as equações de força, o solver balanceado, a
lei ajustada, os limiares ou os artefatos T01–T11.1. A e D são lidos dos CSVs
congelados da T08 e auditados pelas APIs públicas com
`rtol=5e-12, atol=5e-14` antes da primeira solução E. A T12 também não executa
uma nova varredura nem ajusta \(C\) ou \(p\).

Os parâmetros físicos são

\[
a=E_0=1,\qquad ka=0.1,\qquad f_0=0,
\]

com os valores positivos de \(f_1\) e as geometrias congeladas na T08.

## Referência e métricas

A referência principal é a força de interação completa,

\[
\mathbf F^E=\mathbf F^E_{\mathrm{int}}
=\mathbf F^E_{\mathrm{ext-sc}}+\mathbf F^E_{\mathrm{ss}},
\]

e não a força total. A e D recebem explicitamente componente \(z=0\), mas a
componente \(F_z\) calculada pelo Modelo E participa de todas as métricas.

Para vetores tridimensionais,

\[
\mathcal R(\mathbf F)=
\left[\frac1N\sum_i\lVert\mathbf F_i\rVert_2^2\right]^{1/2}.
\]

Os erros principais são

\[
\varepsilon_A^E=
\frac{\mathcal R(\mathbf F^A-\mathbf F^E)}{\mathcal R(\mathbf F^E)},
\qquad
\varepsilon_D^E=
\frac{\mathcal R(\mathbf F^D-\mathbf F^E)}{\mathcal R(\mathbf F^E)}.
\]

Também são registrados os erros contra
\(\mathbf F^E_{\mathrm{ext-sc}}\), o erro simétrico RMS e a identidade vetorial

\[
\mathbf F^E-\mathbf F^A=
(\mathbf F^D-\mathbf F^A)
+(\mathbf F^E_{\mathrm{ext-sc}}-\mathbf F^D)
+\mathbf F^E_{\mathrm{ss}}.
\]

As amplitudes \(X_{D-A}\), \(X_{\mathrm{Mie/ext-sc}}\) e \(X_{\mathrm{ss}}\)
são normalizadas por \(\mathcal R(\mathbf F^E)\). Elas são amplitudes RMS
não negativas, mas não são parcelas escalares aditivas: a reconstrução ocorre
com vetores assinados e admite reforço ou cancelamento.

Razões não aplicáveis são armazenadas como zero finito junto a uma flag e uma
razão textual; nenhum piso dimensional é introduzido.

## Lei e limiares congelados

A previsão usada sem novo ajuste é

\[
\widehat\varepsilon_A=2.6353684041458636\,
\rho_1^{1.1088518115798773}.
\]

Os limiares congelados são:

| tolerância | limiar de \(\rho_1\) |
|---:|---:|
| 0.01 | 0.0053990295322641655 |
| 0.05 | 0.02000077753569526 |
| 0.10 | 0.03914887870730305 |

## Convergência e diagnósticos

Cada sentinela é resolvido com o Modelo E para
\(L_{\max}=2,3,\ldots,13\). O cálculo chega no mínimo a \(L_{\max}=5\) e só
pode parar quando total, interação, external–scattered e
scattered–scattered apresentam, cada um, duas mudanças sucessivas aplicáveis
menores ou iguais a \(10^{-5}\). Um canal sem confirmação até 13 é registrado
como `unconfirmed`, nunca como divergente.

Em cada ordem são preservados \(\kappa_2(A_q)\), erro retroativo balanceado,
fechamentos \(b=a+Ud\) e \(d=Db\), resíduo da decomposição, máximo \(|F_z|\),
dimensão reduzida, forças dos cinco canais e flags de finitude. Uma solução
aceita deve usar `balanced_sqrt`, satisfazer \(\kappa_2(A_q)<10\), manter os
erros físicos abaixo de \(10^{-12}\) e respeitar a simetria em \(z\).

## Elegibilidade, auditoria e gate

`threshold_metric_applicable` exige interação E confirmada e escala resolvida.
`mechanism_decomposition_applicable` exige confirmação da interação e dos dois
canais de mecanismo. `prediction_metric_applicable` exige ainda erro observado
positivo.

A auditoria dos limiares reporta globalmente, por \(N\) e por família:
elegíveis, seguros previstos, seguros observados, falsos seguros, falsos
inseguros, pior erro previsto como seguro e maior \(\rho_1\) observado como
seguro. A lei congelada é avaliada por RMSE logarítmico, fatores mediano,
percentil 90 e máximo, fração dentro de fator 2 e Spearman; nenhuma rotina de
regressão é chamada.

O gate interno para uma futura T13 requer: ao menos 80% de interação E
confirmada; todos os diagnósticos aprovados; RMSE logarítmico no máximo
\(\ln2\); ao menos 80% dentro de fator 2; cobertura prevista no limiar de 5%;
e nenhum falso seguro em 5%. Um gate falso não invalida a conclusão
computacional da T12 e não autoriza alterar o manifesto, tolerâncias ou ajuste.

## Artefatos e aceite

São produzidos deterministicamente quatro CSVs e uma figura 2×2:

- `t12_sentinel_manifest.csv`;
- `t12_model_e_convergence.csv`;
- `t12_model_comparison.csv`;
- `t12_threshold_audit.csv`;
- `t12_model_e_sentinel_audit.png`.

O CSV bruto de convergência é publicado atomicamente somente após a campanha
completa. `--analyze-only` valida sua completude e regenera apenas resultados
derivados. O aceite exige os 28 IDs exatos, auditorias A/D aprovadas,
diagnósticos finitos, testes sem warnings, duas análises byte-idênticas,
preservação byte a byte de T01–T11.1 e inspeção visual da figura.

## Limitações

A conclusão permanece restrita a esferas idênticas e fixas, fluido ideal,
plano nodal, \(ka=0.1\), contrastes positivos amostrados, sete estratos
\(N\leq4\) e força de interação do Modelo E. Explicitamente,

\[
\boxed{
\text{aprovação sentinela em }N\leq4
\ne
\text{validação externa em }N=6,10
\ne
\text{critério universal}.
}
\]

T13 e T14 não são iniciadas por esta tarefa.
