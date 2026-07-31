# T12.3 — critério mecanístico com validação agrupada

## Objetivo

A T12.3 testa, de forma confirmatória e sem novos solves, se um descritor
mecânico mais completo que a lei isolada em \(\rho_1\) pode delimitar regiões
nas quais o erro do Modelo A contra o Modelo E fica abaixo de 1%, 5% ou 10%.
O conjunto é exatamente o dos 28 sentinelas confirmados de T12.1, com
\(N\in\{2,3,4\}\). Os casos externos \(N=6,10\) não são lidos, inspecionados
ou usados.

O commit-base é
`87f67452098c2a85c15c21d26abe9468fa875776`, que encerrou T12.2 com
`NO_GO_T13_RHO1_NOT_QUANTITATIVE`. Nenhuma força, solver, grupo, resposta ou
definição de preditor é alterada.

## Dados e grupos congelados

A fonte exclusiva de respostas é
`results/data/t12_1_resolved_comparison.csv`. Ela contém 28 IDs únicos, em
ordem canônica, distribuídos nos sete grupos:

- `n2_pair`;
- `n3_compact`, `n3_irregular` e `n3_linear`;
- `n4_compact`, `n4_irregular` e `n4_linear`.

Todas as respostas de força de interação do Modelo E estão diretamente
confirmadas. A análise preserva o alvo positivo \(\varepsilon_A^E\),
\(\Lambda_{\max}\), \(\rho_1\), \(f_1\), a separação e os grupos publicados.

O descritor geométrico é

\[
\Lambda_{\max}
=|f_1|\max_i\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3.
\]

Ele é uma soma geométrica máxima de acoplamentos dipolares inverso-cúbicos,
não uma norma, autovalor, valor singular ou raio espectral. Por contraste,
\(\rho_1\) é o raio espectral do operador dipolar balanceado e retém fase e
topologia espectral do reespalhamento. Os dois são mecanisticamente
relacionados, mas não equivalentes.

## Modelos pré-especificados

P0 e P3 são baselines congelados e nunca reajustados:

\[
P0=2.6353684041458636\,\rho_1^{1.1088518115798773},
\]

\[
P3=14.73950709797405\,\rho_1^{1.4226504975598322}.
\]

M1 é testado primeiro:

\[
\log\widehat\varepsilon
=\beta_0+\beta_\Lambda\log\Lambda_{\max},
\qquad
\widehat\varepsilon=C_\Lambda\Lambda_{\max}^{\alpha_\Lambda}.
\]

M2 é apenas o segundo nível pré-especificado:

\[
\log\widehat\varepsilon
=\beta_0+\beta_\Lambda\log\Lambda_{\max}
+\beta_\rho\log\rho_1,
\]

\[
\widehat\varepsilon
=C_{\Lambda\rho}\Lambda_{\max}^{\alpha_\Lambda}
\rho_1^{\alpha_\rho}.
\]

Ambos usam mínimos quadrados ordinários não ponderados em logaritmos naturais,
com intercepto, sem regularização, seleção de variáveis, descarte de casos ou
transformações posteriores.

## LOGO externo e calibração aninhada

Cada um dos sete grupos é retirado integralmente. O modelo é ajustado nos 24
casos dos outros seis grupos e prevê os quatro casos externos. Assim, cada
sentinela recebe exatamente uma previsão OOF por candidato.

Dentro de cada treino externo, outra LOGO nos seis grupos disponíveis produz
24 previsões internas honestas. O fator conservador é

\[
s^{(-g)}=\exp\left[
\max_{j\notin g}
\left(\log\varepsilon_j-\log\widehat\varepsilon_j^{\mathrm{inner}}
\right)
\right].
\]

A previsão de segurança externa é

\[
\widehat\varepsilon_{\mathrm{safe}}
=s^{(-g)}\widehat\varepsilon_{\mathrm{OOF}}.
\]

Um caso só é previsto seguro quando
\(\widehat\varepsilon_{\mathrm{safe}}<\tau\). A igualdade é insegura. A
resposta observada é segura apenas quando \(\varepsilon_A^E<\tau\). P0 e P3
são funções fixas; seus fatores usam somente resíduos dos 24 casos de treino
externo, o equivalente honesto da LOGO interna sem qualquer reajuste.

Os mínimos de antivacuidade são 3, 8 e 12 casos previstos seguros em 1%, 5% e
10%, respectivamente.

## Diagnósticos e gate

As métricas OOF pontuais são RMSE e MAE logarítmicos, frações dentro de fatores
2 e 1,5, Spearman e pior razão multiplicativa. Elas são calculadas globalmente
e por grupo, \(N\), família, \(f_1\) e separação. O bootstrap reamostra grupos
inteiros com semente 1203 e 10.000 amostras; seus intervalos de 95% são apenas
descritivos. A influência leave-one-case-out não remove casos.

M2 é diagnosticado por correlação dos log-preditores, posto, número de condição
da matriz padronizada e estabilidade dos sinais. A etiqueta
`UNSTABLE_COLLINEARITY` é aplicada se houver singularidade, condição acima de
\(10^3\) ou se os dois expoentes mecanicamente esperados positivos aparecerem
em menos de 80% dos folds. Não se usa pseudoinversa customizada ou
regularização.

M1 passa somente com zero falso seguro nas três tolerâncias, antivacuidade,
RMSE OOF menor ou igual a \(\ln2\), pelo menos 85% dentro de fator 2, Spearman
de pelo menos 0,90, todos os folds válidos, expoente de
\(\Lambda_{\max}\) positivo e integridade congelada. M2 só poderia substituir
M1 se M1 falhasse e se, além dos critérios comuns, fosse identificável, tivesse
sinais estáveis e melhorasse o RMSE em pelo menos 5% ou ampliasse estritamente
a cobertura segura sem perdas e sem falsos seguros.

## Resultados

Os ajustes completos, rotulados como descritivos e excluídos do gate, são

\[
\widehat\varepsilon_{M1}
=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\]

\[
\widehat\varepsilon_{M2}
=5.0396777007270535\,
\Lambda_{\max}^{1.2602714475189609}
\rho_1^{0.13234182295409233}.
\]

| Modelo | RMSE log OOF | MAE log | fator 2 | fator 1,5 | Spearman | pior razão |
|---|---:|---:|---:|---:|---:|---:|
| P0 | 0.810603355799503 | 0.633003136612883 | 0.607142857142857 | 0.428571428571429 | 0.970990695128626 | 15.1576399663195 |
| P3 | 0.577854414091078 | 0.336406995083200 | 0.928571428571429 | 0.750000000000000 | 0.970990695128626 | 13.2499440166714 |
| M1 | 0.629389092024730 | 0.400097164805058 | 0.892857142857143 | 0.714285714285714 | 0.946907498631637 | 14.5542179607905 |
| M2 | 0.662457567037336 | 0.426097483531842 | 0.928571428571429 | 0.678571428571429 | 0.943623426382047 | 17.3299357001287 |

As previsões P0 e P3 são avaliações das leis congeladas. Como não há ajuste
novo nelas, sua natureza OOF decorre de não dependerem de nenhum caso da
T12.3; somente a margem conservadora respeita a separação externa.

| Modelo | \(\tau\) | previstos seguros | seguros observados | falsos seguros | falsos inseguros | cobertura segura |
|---|---:|---:|---:|---:|---:|---:|
| P3 | 1% | 7 | 7 | 0 | 0 | 1.000000000000000 |
| P3 | 5% | 13 | 14 | 0 | 1 | 0.928571428571429 |
| P3 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |
| M1 | 1% | 4 | 7 | 0 | 3 | 0.571428571428571 |
| M1 | 5% | 9 | 14 | 0 | 5 | 0.642857142857143 |
| M1 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |
| M2 | 1% | 4 | 7 | 0 | 3 | 0.571428571428571 |
| M2 | 5% | 9 | 14 | 0 | 5 | 0.642857142857143 |
| M2 | 10% | 14 | 20 | 0 | 6 | 0.700000000000000 |

M1 supera todos os mínimos 3/8/12 e não produz falso seguro conservador. O
antigo falso seguro `n2_pair_f0.8_d2.5` tem erro observado
0.12057318984999543; sua previsão M1 pontual é 0.07162592847780698, mas a
margem aninhada do fold eleva a previsão segura para 0.15255694330602437.
Logo ele é corretamente classificado como inseguro em 10%.

Em M2, a correlação completa entre \(\log\Lambda_{\max}\) e
\(\log\rho_1\) é 0.995454787228431 e a condição padronizada é
20.95288536113825. Embora longe de \(10^3\), \(\alpha_\rho\) é positivo em
apenas quatro dos sete folds; o bootstrap também atravessa zero para ambos os
expoentes. M2 recebe `UNSTABLE_COLLINEARITY`, piora o RMSE OOF em 5,25% frente
a M1 e não amplia a cobertura conservadora.

## Decisão

M1 passa literalmente os oito itens do gate. Pela ordem hierárquica, M2 não
substitui M1. A decisão oficial é

```text
GO_T13_VALIDATE_LAMBDA_MAX
```

Isso libera somente uma tarefa posterior de validação externa em \(N=6,10\).
Não constitui validação externa, universalidade ou aprovação para clusters
maiores.

## Produtos e verificação

Foram criados o módulo `mechanistic_validity.py`, o script
`analyze_t12_3_mechanistic_validity.py`, dois arquivos de teste, oito CSVs e
uma figura. Duas execuções integrais no mesmo ambiente produziram bytes
idênticos. O ambiente foi Python 3.12.3, NumPy 2.5.1, SciPy 1.18.0 e
Matplotlib 3.11.1.

Os hashes oficiais dos nove novos artefatos são registrados em
`docs/HANDOFF.md`. Todos os 61 artefatos versionados anteriores — os 54
congelados antes da T12.2 mais os sete produtos da T12.2 — permaneceram
byte a byte inalterados.

## Limitações

O conjunto tem apenas 28 casos internos, \(N\leq4\), \(ka=0.1\),
\(f_0=0\), contrastes positivos amostrados, famílias planares fixas, esferas
idênticas e força completa de interação do Modelo E. Os intervalos bootstrap
são diagnósticos; a correlação forte limita a interpretação de M2. Nenhum
novo solve, caso \(N=6,10\), T13 ou T14 foi executado.
