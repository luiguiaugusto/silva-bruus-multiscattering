# T14.1 — confirmação em grande \(N\) do critério congelado \(\Lambda_{\max}\)

## Objetivo

Esta tarefa testa, sem recalibração, a transferência do critério mecanicista
congelado nas T12.3–T14 para \(N=45\) e \(N=105\). O observável confirmatório é
o erro entre a força pairwise A e a força completa de interação do Modelo E,

\[
\varepsilon_A^E=
\frac{\operatorname{RMS}(\mathbf F_A-\mathbf F_E)}
{\operatorname{RMS}(\mathbf F_E)}.
\]

A seleção geométrica, as previsões e o gate são congelados antes de consultar
qualquer resposta do Modelo E nas coordenadas T14.1. A tarefa não ajusta novo
modelo, não muda M1/P3 e não inicia a T15.

## Base e cronologia cega

O commit-base é
`7c6329ffa6208ac517602f4e411172128fc43465`. A execução possui duas fases:

1. Fase A: criar código completo, geometrias, protocolo, previsões e hashes;
   validar, fazer commit e push com
   `chore: preregister T14.1 large-N confirmation`;
2. Fase B: somente depois desse push, executar o Modelo E, analisar, auditar,
   documentar, fazer commit e push com
   `feat: validate frozen lambda max at large N`.

Nenhum solve do Modelo E T14.1 pode ocorrer na Fase A. É permitido calcular
\(\rho_1\) em \(L_{\max}=1\) para congelar o diagnóstico P3 secundário.

## Regime físico

Mantêm-se \(a=1\), \(E_0=1\), \(ka=0.1\), \(f_0=0\), \(f_1=0.8\), centros
planares, esferas idênticas e não sobrepostas, fluido ideal e o Modelo E
completo já aprovado. Não há mudança em solver, precisão, força ou critério.

## Amostra externa

São 24 casos na ordem literal congelada:

- primeiro \(N=45\), depois \(N=105\);
- em cada tamanho: `linear`, `compact`, `irregular`;
- em cada família: níveis 1, 2, 3 e 4.

As famílias compactas usam malhas triangulares \(q=9\to45\) e
\(q=14\to105\). A família irregular reutiliza exatamente a perturbação
determinística da T14, com amplitude 0.15, recentralização e renormalização.
Cada forma fixa é dilatada uniformemente para um dos quatro alvos congelados:

\[
\Lambda_{\max}\in\{
0.0031111241226691642,
0.011108933664494051,
0.025457132710914911,
0.065350897425260762
\}.
\]

## Acoplamento local congelado

Para cada partícula,

\[
\Lambda_i=|f_1|\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3,
\qquad
\Lambda_{\max}=\max_i\Lambda_i.
\]

O pré-registro guarda o vetor completo em ordem, SHA-256 e as estatísticas:
mínimo, média, mediana, desvio padrão, percentis 10 e 90, máximo,
\(\overline\Lambda/\Lambda_{\max}\), fração com
\(\Lambda_i\ge0.9\Lambda_{\max}\) e primeiro índice do máximo.

## Predições congeladas

M1 usa \(\Lambda_{\max}\); P3 usa \(\rho_1\) apenas como diagnóstico. Os
coeficientes, expoentes, fatores conservadores e limiares são lidos das APIs
congeladas. As contagens cegas obrigatórias de M1 em 1%, 5% e 10% são
6/12/18. M2 permanece excluído por colinearidade instável. P3 não pode mudar a
decisão confirmatória de M1.

## Protocolo Modelo E

Para cada caso, calcular \(L_{\max}=2,3,\ldots,13\), sem usar apenas ordens
ímpares e sem parar antes de \(L_{\max}=5\). Total, interação,
external–scattered e scattered–scattered exigem duas mudanças RMS sucessivas,
aplicáveis e menores ou iguais a \(10^{-5}\). Só há parada antecipada quando os
quatro canais são confirmados.

O runner usa por padrão um worker e uma thread BLAS, executa \(N=45\) antes de
\(N=105\), faz estimativa conservadora de memória antes de cada ordem e salva
cache atômico retomável somente em `/tmp`. O CSV bruto oficial é publicado
atomicamente após a campanha completa ou a consolidação explícita de limite de
recurso. Não se executa ordem acima de 13.

Cada ordem registra dimensões, memória prevista e observada, tempos,
condicionamento balanceado, erro retroativo, fechamentos, resíduo da
decomposição, finitude, planaridade e mudanças dos quatro canais. O gate
numérico exige solver `balanced_sqrt`, condicionamento menor que 10 e todos os
erros/fechamentos menores que \(10^{-12}\).

## Elegibilidade

Um caso é elegível somente quando a campanha terminou, a interação foi
confirmada, o registro final passou todos os diagnósticos, o erro A–E é
aplicável e toda a identidade da Fase A foi preservada. Não há imputação,
extrapolação nem piso absoluto. Casos inelegíveis são mantidos com motivo.

## Análises

Para M1 e P3, reportam-se métricas multiplicativas globais, por tamanho, por
família e por nível. A auditoria conservadora registra seguros previstos e
observados, falsos seguros, falsos inseguros e IDs em 1%, 5% e 10%.

A transferência pareada calcula

\[
R_{105/45}=\frac{\varepsilon_A^E(105)}{\varepsilon_A^E(45)},
\quad
R_{45/28},
\quad
R_{105/28},
\]

e combina, apenas de forma descritiva, a sequência
\(N=15\to28\to45\to105\). Também se relaciona o erro com a estrutura de
\(\Lambda_i\), o diâmetro acústico e o custo computacional, sem substituir M1.

## Diagnóstico de tendência em grande N

Com pelo menos dez pares aplicáveis:

- não há deterioração sistemática se a mediana de \(R_{105/45}\) for
  \(\le1.10\) e o percentil 90 for \(\le1.25\);
- há deterioração sistemática se a mediana for \(>1.25\) e pelo menos 9 de 12
  razões forem \(>1.10\);
- os demais casos são classificados como tendência mista;
- com menos de dez pares, o diagnóstico é inconclusivo.

Esse diagnóstico não altera o gate confirmatório.

## Gate congelado

A suficiência exige pelo menos 20/24 elegíveis, 10 por tamanho, 6 por família,
5 por nível, todos os 6/12/18 casos cegamente seguros elegíveis, integridade dos
cinco artefatos da Fase A, integridade de todos os artefatos anteriores e
nenhuma ordem acima de 13.

Falha de convergência/diagnóstico produz:

```text
INCONCLUSIVE_T14_1_INSUFFICIENT_MODEL_E_CONVERGENCE
HOLD_T15_T14_1_INCONCLUSIVE
```

Limitação computacional produz:

```text
INCONCLUSIVE_T14_1_RESOURCE_LIMIT
HOLD_T15_T14_1_INCONCLUSIVE
```

Com suficiência completa, M1 passa apenas com zero falso seguro nas três
tolerâncias, contagens 6/12/18, RMSE log global e por tamanho \(\le\ln2\),
fração em fator 2 global \(\ge0.80\) e por tamanho \(\ge0.75\), Spearman global
\(\ge0.90\) e protocolo imutável. O resultado é um dos literais:

```text
PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX
GO_T15_SYNTHESIS_AND_MANUSCRIPT
```

ou

```text
FAIL_T14_1_LARGE_N_FROZEN_LAMBDA_MAX
GO_T15_SYNTHESIS_WITH_LARGE_N_BREAKDOWN
```

## Artefatos

A Fase A publica manifesto, acoplamento local, previsões, protocolo e hashes
anteriores. A Fase B publica CSV bruto de convergência, forças, resumo de casos,
predições reveladas, métricas, auditoria de limiares, pares, sequência combinada,
desempenho, gate e figura de seis painéis.

## Auditoria e aceite

Após a revelação, `--audit-existing` recalcula seis casos estratificados nas
ordens finais e compara forças, canais, \(\Lambda_i\), \(\Lambda_{\max}\),
\(\rho_1\) e diagnósticos. Duas análises consecutivas devem ser byte-idênticas.
Todos os testes devem passar com warnings como erros, os 95 artefatos anteriores
devem permanecer byte a byte, documentos não podem conter caracteres ASCII de
controle e a figura deve ser inspecionada visualmente.

## Limitações

As conclusões se restringem a duas realizações determinísticas grandes, três
famílias planares, quatro níveis de acoplamento, \(ka=0.1\), \(f_1=0.8\) e ao
Modelo E aprovado. Não demonstram universalidade, não recalibram M1 e não
autorizam extrapolação além da amostra.
