# Prompt da T08 — Transferibilidade para N ≤ 10 e critério de acoplamento coletivo

Execute a T08 do projeto `silva-bruus-multiscattering`.

# T08 — Transferibilidade para N <= 10 e critério de acoplamento coletivo

Esta é a última tarefa computacional planejada para o artigo. Seu objetivo não é ampliar indefinidamente o número de partículas, mas testar se os resultados obtidos para N <= 4 permitem prever o erro da aproximação pairwise em clusters com N = 6 e N = 10.

A T08 deve terminar com o congelamento dos dados do artigo, independentemente de o critério coletivo funcionar ou falhar.

---

## 1. Estado inicial obrigatório

Trabalhe sobre a `main` remota no commit:

`bd3752ef2d4fa61ec10887445575f507d8c8cd6b`

Mensagem desse commit:

`feat: add multipolar Model D convergence study`

Antes de alterar qualquer arquivo:

1. leia integralmente:
   - `AGENTS.md`;
   - `README.md`;
   - `TASKS.md`;
   - `docs/CONVENTIONS.md`;
   - `docs/DECISIONS.md`;
   - `docs/HANDOFF.md`;
   - `TAREFA_T07_MODELO_D_CONVERGENCIA_MULTIPOLAR.md`;

2. execute:
   - `git status --short`;
   - `git branch --show-current`;
   - `git rev-parse HEAD`;
   - `git rev-parse origin/main`;

3. confirme:
   - branch `main`;
   - árvore limpa;
   - `HEAD == origin/main`;
   - commit-base exatamente igual ao hash acima.

Não descarte nem sobrescreva alterações preexistentes. Se a árvore estiver suja ou a `main` estiver em outro commit, pare e informe.

Antes das modificações, registre SHA-256 de todos os artefatos já existentes em `results/`. Todos deverão permanecer byte a byte inalterados.

---

## 2. Pergunta científica

A T08 deve responder:

> Um parâmetro escalar de acoplamento coletivo, calibrado apenas com clusters de N <= 4, consegue prever a validade da aproximação pairwise em clusters não usados na calibração, com N = 6 e N = 10?

O erro primário é o erro de Silva–Bruus em relação ao Modelo D multipolarmente convergido.

A análise secundária deve medir o erro residual depois de todas as interações de dois corpos terem sido calculadas no mesmo truncamento multipolar da solução global.

Não ajuste um critério usando N = 6 ou N = 10 e depois chame isso de validação. Esses tamanhos constituem o conjunto de teste externo.

---

## 3. Escopo físico fixo

Use exclusivamente:

- plano nodal;
- partículas esféricas, idênticas e fixas;
- fluido ideal, homogêneo, ilimitado e invíscido;
- `ka = 0.1`;
- `radius = 1`;
- `k = 0.1`;
- `energy_density = 1`;
- `f0 = 0`;
- `f1 in {0.1, 0.4, 0.8, 1.0}`;
- `d_min/a in {2.1, 2.5, 3.0, 4.0, 6.0, 10.0}`;
- `N in {2, 3, 4, 6, 10}`;
- força external–scattered já validada;
- coeficientes Rayleigh dominantes de cada multipolo, exatamente como na T07.

Não incluir:

- `N > 10`;
- partículas fora do plano nodal;
- partículas diferentes;
- `ka` adicional;
- contrastes negativos;
- ensembles aleatórios;
- trajetórias ou dinâmica;
- forças scattered–scattered;
- T-matrix exata de Mie;
- viscosidade;
- streaming;
- paredes;
- contato;
- torque;
- FMM;
- decomposição conectada acima de quatro corpos.

Não executar novamente as varreduras completas de 1.920 casos da T05 ou da T06.

---

## 4. Modelos usados na T08

### 4.1 Modelo A

Preserve exatamente o Modelo A existente:

\[
\mathbf F_i^A
=
\sum_{j\ne i}\mathbf F_{ij}^{\mathrm{SB}}.
\]

Não altere a implementação ou a definição de Silva–Bruus.

### 4.2 Baseline pairwise multipolarmente compatível

Defina uma nova API, sem alterar o significado histórico do Modelo B, para

\[
\mathbf F_i^{B_L}
=
\sum_{j\ne i}
\mathbf F_{ij}^{D,N=2,L}.
\]

Para cada par não ordenado:

1. resolva o dímero isolado com o Modelo D;
2. use o mesmo \(L\) empregado na solução global;
3. embuta as duas forças nas posições correspondentes do cluster;
4. some todos os pares.

Essa grandeza deve ser chamada na documentação de:

`multipolarly matched pairwise baseline`

ou

`matched pairwise baseline B_L`.

Ela não constitui um novo solver físico. É apenas a soma pairwise de soluções de dímeros do Modelo D.

Propriedades obrigatórias:

\[
B_{L=1}=B
\]

para o Modelo B histórico, e

\[
B_L=D_L
\]

para \(N=2\), dentro da precisão numérica.

Use cache determinístico por separação, `f1` e `L` para não resolver repetidamente dímeros equivalentes.

### 4.3 Modelo D

Use sem alterações:

`solve_multipolar_nodal_interaction_forces`.

A solução global será denotada por

\[
\mathbf F^{D_L}.
\]

A identidade vetorial

\[
\mathbf F^{D_L}-\mathbf F^A
=
\left(\mathbf F^{B_L}-\mathbf F^A\right)
+
\left(\mathbf F^{D_L}-\mathbf F^{B_L}\right)
\]

deve ser verificada programaticamente.

Não use o Modelo C como referência final. O resultado em \(L=1\) deve ser mantido apenas como diagnóstico da correção multipolar e como regressão com as etapas anteriores.

---

## 5. Geometrias

Todas as geometrias devem:

- estar em `z = 0`;
- ser centradas no centroide;
- possuir distância mínima exatamente igual ao parâmetro `d_min`;
- ser determinísticas;
- não usar números aleatórios.

Crie um novo módulo para as famílias da T08. Não modifique as geometrias científicas anteriores.

### 5.1 Par, N = 2

Use apenas uma família denominada `pair`.

### 5.2 Cadeia linear

Para qualquer \(N\):

\[
x_i =
\left(i-\frac{N-1}{2}\right)d_{\min},
\qquad
y_i=z_i=0.
\]

Nome da família: `linear`.

### 5.3 Cluster compacto

Use:

- N = 3: `equilateral_trimer`;
- N = 4: `square_quartet`;
- N = 6: patch triangular com linhas 3–2–1;
- N = 10: patch triangular com linhas 4–3–2–1.

Para N = 6, use o template bidimensional:

\[
(0,0),(1,0),(2,0),
(1/2,\sqrt3/2),(3/2,\sqrt3/2),
(1,\sqrt3).
\]

Para N = 10:

\[
\begin{aligned}
&(0,0),(1,0),(2,0),(3,0),\\
&(1/2,\sqrt3/2),(3/2,\sqrt3/2),(5/2,\sqrt3/2),\\
&(1,\sqrt3),(2,\sqrt3),\\
&(3/2,3\sqrt3/2).
\end{aligned}
\]

Centre e escale cada template para `d_min`.

Nome da família: `compact`.

### 5.4 Cluster irregular

Use:

- N = 3: `scalene_trimer`;
- N = 4: `irregular_quartet`.

Para N = 6, use o template bruto:

```python
[
    [0.00, 0.00],
    [1.05, 0.08],
    [2.18, -0.06],
    [0.32, 1.14],
    [1.49, 1.03],
    [0.91, 2.22],
]
```

Para N = 10:

```python
[
    [0.00, 0.00],
    [1.10, 0.05],
    [2.25, -0.10],
    [3.30, 0.18],
    [0.35, 1.12],
    [1.48, 0.95],
    [2.62, 1.30],
    [0.82, 2.15],
    [1.92, 2.30],
    [1.28, 3.28],
]
```

Para cada template:

1. calcule a menor distância entre os pontos;
2. normalize essa distância para 1;
3. subtraia o centroide;
4. multiplique por `d_min`;
5. acrescente `z = 0`.

Nome da família: `irregular`.

### 5.5 Número exato de configurações

A enumeração geométrica deve conter:

- 1 configuração para N = 2;
- 3 para N = 3;
- 3 para N = 4;
- 3 para N = 6;
- 3 para N = 10.

Logo:

\[
13\times4\times6=312
\]

configurações físicas.

Não adicione outras geometrias.

---

## 6. Protocolo de convergência multipolar

Para cada configuração, avalie sucessivamente:

\[
L\in\{1,3,5,7,9,11\}.
\]

Somente para \(N\le4\), permita \(L=13\) quando a convergência ainda não estiver confirmada em \(L=11\).

Para \(N=6\) e \(N=10\), o limite absoluto é \(L=11\). Casos não convergidos até esse ponto devem ser registrados como não confirmados. Não aumente silenciosamente para \(L=13,15,\ldots\).

Em cada ordem, calcule:

\[
\mathbf R_L
=
\mathbf F^{D_L}-\mathbf F^{B_L}.
\]

Para duas ordens sucessivas, defina:

\[
\delta_D(L)
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^{D_L}-\mathbf F^{D_{L-2}}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^{D_L}\right)},
\]

\[
\delta_B(L)
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^{B_L}-\mathbf F^{B_{L-2}}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^{D_L}\right)},
\]

\[
\delta_R(L)
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf R_L-\mathbf R_{L-2}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^{D_L}\right)}.
\]

Uma quantidade estará convergida somente quando suas duas últimas diferenças sucessivas forem menores ou iguais a

\[
10^{-3}.
\]

Registre separadamente:

- `total_converged`;
- `matched_pairwise_converged`;
- `collective_residual_converged`;
- `joint_converged`;
- `reference_lmax`;
- `maximum_allowed_lmax`;
- últimas duas diferenças de cada quantidade.

Não confunda convergência da força total com convergência do residual coletivo.

Para a análise de \(\varepsilon_A\), exija pelo menos `total_converged`.

Para a análise de \(\varepsilon_B\), exija `joint_converged`.

Defina ainda:

\[
u_R=\max\left[\delta_R(L),\delta_R(L-2)\right].
\]

O residual coletivo só será considerado numericamente resolvido quando

\[
\varepsilon_B>5u_R.
\]

Registre isso como `collective_residual_resolved`.

Casos não convergidos podem ter seus últimos valores disponíveis registrados, mas:

- não podem entrar nos ajustes;
- não podem definir limiares;
- devem aparecer nas figuras como pontos abertos ou cruzes;
- devem ser chamados de `unconfirmed`, nunca de divergentes.

Implemente tratamento scale-aware caso a força de referência seja numericamente nula. Não grave `NaN` ou `inf`; use campos de aplicabilidade explícitos.

---

## 7. Métricas principais

Use:

\[
\varepsilon_A
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^A-\mathbf F^{D}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^{D}\right)},
\]

\[
\varepsilon_B
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^{B_L}-\mathbf F^{D}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^{D}\right)}.
\]

Calcule também:

\[
Y_{\mathrm{2B}}
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^{B_L}-\mathbf F^A
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^D\right)},
\]

\[
Y_{\mathrm{coll}}
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^D-\mathbf F^{B_L}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^D\right)},
\]

e a correção multipolar total:

\[
Y_{\mathrm{mp}}
=
\frac{
F_{\mathrm{RMS}}\!\left(
\mathbf F^{D}-\mathbf F^{D_{L=1}}
\right)}
{F_{\mathrm{RMS}}\!\left(\mathbf F^D\right)}.
\]

Os CSVs devem armazenar essas grandezas como frações, não porcentagens.

---

## 8. Preditores candidatos

Compare exatamente três preditores.

### 8.1 Acoplamento local mínimo

\[
\eta
=
|f_1|
\left(\frac{a}{d_{\min}}\right)^3.
\]

### 8.2 Máximo acoplamento geométrico coletivo

Reutilize a função existente:

\[
\Lambda_{\max}
=
|f_1|
\max_i
\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3.
\]

Não altere sua implementação anterior.

### 8.3 Raio espectral do operador de reespalhamento em L = 1

Use a matriz balanceada da solução em \(L=1\):

\[
\mathbf A_b^{(1)}
=
\mathbf I-\mathbf K_b^{(1)}.
\]

Defina:

\[
\rho_1
=
\max_\nu
\left|
\lambda_\nu\left(\mathbf K_b^{(1)}\right)
\right|,
\qquad
\mathbf K_b^{(1)}
=
\mathbf I-\mathbf A_b^{(1)}.
\]

Não use:

- o raio espectral da matriz do sistema;
- o número de condicionamento;
- a matriz bruta mal escalada.

Documente que \(\rho_1\) mede o reespalhamento permitido pela base dipolar \(L=1\), enquanto \(\Lambda_{\max}\) é um estimador geométrico simples.

---

## 9. Calibração e teste de transferibilidade

### 9.1 Conjunto de calibração

Use somente:

\[
N\in\{2,3,4\}.
\]

### 9.2 Conjunto de teste externo

Use somente:

\[
N\in\{6,10\}.
\]

Nenhum dado de N = 6 ou N = 10 pode:

- selecionar o melhor preditor;
- ajustar expoentes;
- ajustar prefatores;
- definir limiares de validade.

### 9.3 Ajustes descritivos

Para cada preditor \(P\), ajuste no conjunto de calibração:

\[
\varepsilon=C P^p
\]

em espaço logarítmico, separadamente para:

- \(\varepsilon_A\);
- \(\varepsilon_B\).

Para \(\varepsilon_B\):

- exclua N = 2, pois \(B_L=D_L\);
- exclua resíduos não resolvidos;
- exclua casos não convergidos.

Registre:

- número de pontos;
- prefator;
- expoente;
- \(R^2\) logarítmico;
- RMSE logarítmico;
- máximo resíduo logarítmico;
- correlação de Spearman.

Os ajustes são empíricos e descritivos, não provas analíticas.

### 9.4 Seleção do preditor sem vazamento

Selecione o preditor primário usando somente \(\varepsilon_A\) no conjunto de calibração.

Faça validação cruzada deixando de fora, uma por vez, cada combinação `(N, family)` do conjunto N <= 4.

O critério de seleção é o menor RMSE logarítmico agregado de validação cruzada.

Em empate numérico dentro de `1e-12`, use a ordem:

1. `lambda_max`;
2. `rho_l1`;
3. `eta`.

Somente depois dessa seleção avalie o preditor escolhido em N = 6 e N = 10.

No holdout, registre:

- RMSE logarítmico;
- mediana do fator multiplicativo de erro;
- percentil 90 desse fator;
- máximo fator;
- fração das previsões dentro de um fator 2;
- correlação de Spearman;
- resultados separados para N = 6 e N = 10;
- resultados agrupados.

Também registre o desempenho dos três preditores para transparência, mas não escolha o vencedor usando o holdout.

---

## 10. Limiares empíricos de validade

Para:

\[
\tau\in\{0.01,0.05,0.10\},
\]

construa, apenas com a calibração, um limiar conservador \(P_\tau\).

Defina \(P_\tau\) como o maior valor observado do preditor para o qual todos os casos de calibração com \(P\le P_\tau\) satisfazem:

\[
\varepsilon_A\le\tau.
\]

Exija pelo menos oito casos de calibração abaixo do limiar. Caso contrário, registre o limiar como não disponível por meio de uma flag, nunca `NaN`.

Aplique o limiar congelado ao holdout e reporte:

- número de casos elegíveis;
- número previsto como seguro;
- cobertura;
- falsos seguros, isto é, casos com \(P\le P_\tau\) e \(\varepsilon_A>\tau\);
- pior erro observado entre os casos previstos como seguros.

Esses limiares devem ser chamados de:

`empirical nodal-plane thresholds within the sampled domain`.

Nunca os apresente como universais.

---

## 11. Critério diagnóstico de sucesso científico

A conclusão física não deve fazer parte dos testes automatizados.

Registre `criterion_supported = true` somente se, para o preditor selecionado:

1. pelo menos 80% do holdout tiver convergência total confirmada;
2. o RMSE logarítmico do holdout for menor ou igual a \(\ln 2\);
3. pelo menos 80% das previsões estiverem dentro de um fator 2;
4. o limiar de 5% tiver pelo menos um caso holdout previsto como seguro;
5. não houver falso seguro no limiar de 5%.

Se essas condições não forem satisfeitas:

- a T08 ainda estará computacionalmente concluída;
- registre `criterion_supported = false`;
- conclua que nenhum critério escalar suficientemente transferível foi demonstrado no domínio testado;
- não amplie N, L, geometrias ou parâmetros.

---

## 12. Implementação

Prefira criar:

```text
src/acoustic_ms/cluster_families.py
src/acoustic_ms/transferability.py
scripts/run_t08_transferability.py
scripts/analyze_t08_transferability.py
tests/test_t08_cluster_families.py
tests/test_t08_transferability.py
tests/test_t08_analysis.py
TAREFA_T08_TRANSFERIBILIDADE_CRITERIO_ACOPLAMENTO.md
```

Atualize `src/acoustic_ms/__init__.py` somente para exportar as novas APIs públicas.

O gerador caro e a análise devem ser separados:

### `run_t08_transferability.py`

Responsável por:

- enumerar exatamente 312 configurações;
- executar a convergência adaptativa;
- calcular A, \(B_L\), D, métricas e preditores;
- gerar os dados brutos;
- manter ordenação determinística;
- oferecer um modo `--audit-existing` que regenere apenas uma amostra estratificada e compare com os CSVs já existentes.

A amostra de auditoria deve incluir pelo menos:

- N = 2, 4 e 10;
- contraste fraco e forte;
- separação próxima e distante;
- as três famílias para N = 10.

### `analyze_t08_transferability.py`

Responsável por:

- ler exclusivamente os CSVs brutos da T08;
- aplicar as regras de elegibilidade;
- fazer calibração, validação cruzada e holdout;
- selecionar o preditor;
- calcular limiares;
- gerar tabelas e figuras.

Não recalcule as 312 configurações durante a análise.

Não altere os solvers científicos da T07 para acelerar a T08. Otimizações devem ficar na nova camada de cache e orquestração.

Se o piloto indicar tempo total excessivo, otimize reutilização de pares e separações antes da rodada principal. Não reduza o grid nem altere os limites de L silenciosamente.

---

## 13. Artefatos obrigatórios

Gere:

```text
results/data/t08_cases.csv
results/data/t08_forces.csv
results/data/t08_convergence.csv
results/data/t08_predictor_fits.csv
results/data/t08_validity_thresholds.csv
results/figures/t08_predictor_comparison.png
results/figures/t08_transferability.png
```

### `t08_cases.csv`

Uma linha por configuração física, totalizando exatamente 312 linhas.

Deve conter pelo menos:

- `case_id`;
- `split`;
- `particle_count`;
- `family`;
- parâmetros físicos;
- número de pares;
- `reference_lmax`;
- flags de convergência;
- flag de resolução coletiva;
- \(\eta\);
- \(\Lambda_{\max}\);
- \(\rho_1\);
- RMS de A, \(B_L\), D;
- \(\varepsilon_A\);
- \(\varepsilon_B\);
- \(Y_{\mathrm{2B}}\);
- \(Y_{\mathrm{coll}}\);
- \(Y_{\mathrm{mp}}\);
- incertezas finais de truncamento;
- resíduo físico máximo;
- condicionamento balanceado máximo;
- condicionamento bruto máximo;
- aplicabilidade das métricas.

### `t08_forces.csv`

Formato longo, com uma linha por partícula e configuração:

- `case_id`;
- índice da partícula;
- posição;
- componentes de A;
- componentes de \(B_L\);
- componentes de D;
- componentes de \(D-B_L\).

### `t08_convergence.csv`

Uma linha por configuração e ordem calculada, contendo as forças RMS, resíduos, condicionamentos e diferenças sucessivas.

### Figuras

`t08_predictor_comparison.png`:

- grade 2 × 3;
- colunas: \(\eta\), \(\Lambda_{\max}\), \(\rho_1\);
- linha superior: \(\varepsilon_A\);
- linha inferior: \(\varepsilon_B\);
- eixos logarítmicos;
- calibração e holdout visualmente distintos;
- N e famílias distinguíveis;
- casos não confirmados como símbolos abertos ou cruzes;
- ajuste somente da calibração;
- nenhum caso não elegível deve alimentar a linha de ajuste.

`t08_transferability.png`:

- usar somente o preditor selecionado pela validação cruzada;
- mostrar calibração e N = 6/10;
- incluir tolerâncias de 1%, 5% e 10%;
- mostrar ajuste de calibração;
- incluir um painel observado × previsto com diagonal e faixas de fator 2;
- identificar explicitamente que o holdout não participou do ajuste.

As figuras devem ser adequadas para inspeção científica e futura adaptação ao paper.

---

## 14. Testes obrigatórios

Inclua testes para:

1. número correto de partículas;
2. centroide nulo;
3. `z = 0`;
4. distância mínima correta;
5. não sobreposição;
6. determinismo das geometrias;
7. covariância sob rotação e permutação;
8. \(B_{L=1}=B\) histórico;
9. \(B_L=D_L\) para N = 2;
10. identidade vetorial A–\(B_L\)–D;
11. invariância dos preditores sob translação, rotação e permutação;
12. escala \(\Lambda_{\max}\propto d^{-3}\);
13. escala \(\rho_1\) consistente com o operador extraído;
14. raio espectral calculado de \(\mathbf K_b\), não da matriz do sistema;
15. protocolo sintético de duas diferenças sucessivas;
16. exclusão automática de casos não convergidos;
17. exclusão de N = 6/10 da calibração e seleção;
18. exclusão de N = 2 do ajuste de \(\varepsilon_B\);
19. contagem exata de 312 configurações;
20. ordenação determinística;
21. seleção por validação cruzada sem consultar o holdout;
22. limiares calculados apenas na calibração;
23. ausência de `NaN` e `inf`;
24. regressão com a T07.

Para a regressão com a T07:

- leia `t07_cluster_convergence.csv`;
- compare os seis casos canônicos em `f1=0.8`, `d_min/a=2.1`;
- compare as forças do Modelo D nas ordens sobrepostas;
- não altere o CSV da T07;
- use tolerância compatível com ponto flutuante, não igualdade textual entre ambientes diferentes.

Os testes não devem exigir que a hipótese física seja verdadeira. Não faça um teste falhar porque um preditor teve colapso ruim ou porque `criterion_supported` resultou falso.

---

## 15. Validação numérica

Exija:

- todos os resíduos físicos finitos;
- resíduo físico máximo menor que `1e-11`;
- condicionamento balanceado finito;
- nenhum julgamento de convergência baseado no condicionamento bruto;
- forças finitas;
- nenhuma configuração duplicada;
- nenhuma configuração omitida;
- nenhuma sobreposição;
- nenhum caso não convergido incluído em ajuste ou limiar.

Casos próximos, fortes e não convergidos devem ser listados individualmente no relatório final.

Não extrapole resultados de casos não convergidos.

---

## 16. Determinismo sem repetir toda a rodada cara

Execute a varredura completa uma vez.

Depois:

1. registre hashes dos CSVs brutos;
2. execute `--audit-existing` sobre a amostra estratificada;
3. confirme concordância numérica com os CSVs;
4. execute o script de análise duas vezes;
5. confirme que os CSVs derivados e as figuras possuem hashes idênticos nas duas execuções, no mesmo ambiente.

Não execute novamente toda a varredura de 312 casos apenas para testar determinismo.

Registre no `HANDOFF.md`:

- versão do Python;
- NumPy;
- SciPy;
- Matplotlib;
- sistema operacional;
- política de determinismo adotada.

Não espere igualdade byte a byte entre versões diferentes das bibliotecas.

---

## 17. Proteção do trabalho anterior

Nenhum artefato T01–T07 pode ser alterado.

Não modifique os módulos científicos existentes:

```text
src/acoustic_ms/silva_bruus.py
src/acoustic_ms/corrected_pair.py
src/acoustic_ms/force.py
src/acoustic_ms/solver.py
src/acoustic_ms/comparison.py
src/acoustic_ms/scaling.py
src/acoustic_ms/metrics.py
src/acoustic_ms/geometries.py
src/acoustic_ms/multipolar_scattering.py
src/acoustic_ms/multipolar_solver.py
src/acoustic_ms/multipolar_expansion.py
src/acoustic_ms/model_d.py
src/acoustic_ms/translation.py
```

A única exceção entre arquivos existentes do pacote é `src/acoustic_ms/__init__.py`, apenas para exportações.

Atualizações documentais autorizadas:

- `README.md`;
- `TASKS.md`;
- `docs/CONVENTIONS.md`;
- `docs/DECISIONS.md`;
- `docs/HANDOFF.md`.

Não reescreva resultados anteriores como se tivessem sido obtidos na T08.

Corrija no README afirmações antigas de que Modelo D ou N > 4 ainda não existem, preservando o histórico das etapas.

---

## 18. Documentação científica

Documente:

- objetivo da T08;
- definição de \(B_L\);
- diferença entre o Modelo B histórico e o baseline multipolarmente compatível;
- definições de \(\varepsilon_A\), \(\varepsilon_B\), \(Y_{\mathrm{2B}}\), \(Y_{\mathrm{coll}}\) e \(Y_{\mathrm{mp}}\);
- definições de \(\eta\), \(\Lambda_{\max}\) e \(\rho_1\);
- protocolo de convergência;
- separação calibração/holdout;
- prevenção de vazamento;
- validação cruzada;
- definição dos limiares empíricos;
- cobertura de convergência;
- casos não confirmados;
- desempenho no holdout;
- se `criterion_supported` foi verdadeiro ou falso;
- decisão de congelamento dos dados.

Declare explicitamente que a T08 encerra as varreduras do artigo.

As limitações finais devem incluir:

- plano nodal;
- `ka = 0.1`;
- contrastes positivos;
- esferas idênticas;
- partículas fixas;
- N <= 10;
- somente três famílias para N = 6 e 10;
- força external–scattered;
- coeficientes multipolares Rayleigh dominantes;
- ausência de scattered–scattered;
- ausência de viscosidade, streaming, paredes e dinâmica;
- limiares empíricos restritos ao domínio amostrado.

---

## 19. Figuras e inspeção visual

Inspecione visualmente as duas figuras e confirme:

- textos legíveis;
- eixos e títulos não cortados;
- escalas logarítmicas corretas;
- calibração e holdout distinguíveis;
- N e famílias distinguíveis;
- ajustes visíveis;
- tolerâncias de 1%, 5% e 10% visíveis quando aplicáveis;
- diagonal e faixa de fator 2 visíveis;
- casos não confirmados identificáveis;
- nenhuma legenda cobrindo os dados de forma destrutiva;
- ausência de `NaN` ou `inf`;
- ausência de painéis vazios;
- nenhuma afirmação visual de que pontos não convergidos participaram dos ajustes.

---

## 20. Verificações finais

Execute:

```bash
python -m pip install -e ".[dev,plot]"
pytest -q
pytest -q -W error
python scripts/run_t08_transferability.py
python scripts/run_t08_transferability.py --audit-existing
python scripts/analyze_t08_transferability.py
python scripts/analyze_t08_transferability.py
git diff --check
git status --short
git diff --stat
git diff --name-only
```

Faça verificação programática de caracteres ASCII de controle em todos os arquivos alterados, permitindo somente newline e retorno de carro. Não permita tabulações.

Confirme:

- 312 configurações;
- todos os dados finitos;
- hashes anteriores preservados;
- nenhum arquivo protegido alterado;
- nenhum artefato anterior regravado;
- análise determinística;
- figuras inspecionadas;
- nenhuma varredura T05/T06 refeita;
- nenhum N > 10;
- nenhuma decomposição \(\Phi^{(5)}\), \(\Phi^{(6)}\), etc.

---

## 21. Critérios de aprovação

A T08 estará concluída se:

1. todos os testes anteriores continuarem passando;
2. os novos testes passarem com warnings tratados como erros;
3. as 312 configurações forem avaliadas até convergência ou até o limite explícito de L;
4. casos não confirmados forem corretamente marcados e excluídos;
5. \(B_1\) reproduzir o Modelo B;
6. \(B_L=D_L\) para o dímero;
7. a identidade A–\(B_L\)–D fechar numericamente;
8. os três preditores forem calculados;
9. a seleção utilizar somente N <= 4;
10. N = 6/10 forem usados apenas como holdout;
11. os limiares forem construídos apenas na calibração;
12. os resultados físicos, positivos ou negativos, forem documentados sem manipulação;
13. os sete artefatos forem gerados;
14. a auditoria estratificada reproduzir os dados;
15. duas análises sucessivas forem byte-idênticas;
16. os artefatos anteriores permanecerem byte a byte inalterados;
17. nenhum arquivo protegido for alterado;
18. não houver `NaN`, `inf`, tabulações ou caracteres de controle;
19. as figuras forem aprovadas visualmente;
20. a documentação declarar o congelamento do conjunto de dados.

O sucesso da hipótese `criterion_supported` não é condição para conclusão computacional da T08.

---

## 22. Commit e push

Antes do commit, confirme que somente os arquivos autorizados foram modificados.

Use a mensagem:

`feat: test collective-coupling transferability`

Faça commit e push para a `main`.

Não faça push de caches, arquivos temporários ou resultados parciais.

---

## 23. Relatório final

Ao concluir, informe somente:

1. hash completo do commit;
2. confirmação do push;
3. saída de `pytest -q -W error`;
4. arquivos criados e alterados;
5. `git diff --stat`;
6. número exato de configurações;
7. definição implementada de \(B_L\);
8. confirmação \(B_1=B\);
9. confirmação \(B_L=D_L\) para N = 2;
10. erro máximo da identidade A–\(B_L\)–D;
11. fórmulas dos três preditores;
12. tabela de convergência por N e família;
13. lista de casos não confirmados;
14. preditor selecionado pela validação cruzada;
15. tabela dos ajustes de calibração;
16. desempenho no holdout total, N = 6 e N = 10;
17. tabela dos limiares de 1%, 5% e 10%;
18. valor de `criterion_supported`;
19. cobertura de convergência do holdout;
20. resíduos físicos máximos;
21. condicionamentos balanceados máximos;
22. regressão com a T07;
23. hashes dos novos artefatos;
24. confirmação dos hashes anteriores;
25. confirmação da auditoria estratificada;
26. confirmação do determinismo da análise;
27. confirmação da inspeção visual;
28. confirmação da ausência de caracteres de controle;
29. confirmação de que nenhum arquivo protegido foi alterado;
30. confirmação de que nenhuma varredura T05/T06 foi refeita;
31. limitações científicas restantes;
32. confirmação de que os dados computacionais do artigo estão congelados.
