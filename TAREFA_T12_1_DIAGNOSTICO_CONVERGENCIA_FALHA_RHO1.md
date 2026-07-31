# Prompt para execução da T12.1 — Fechamento de convergência e diagnóstico da falha de \(\rho_1\)

Você está trabalhando no repositório:

```text
silva-bruus-multiscattering
```

Execute a **T12.1 — fechamento da convergência multipolar dos sentinelas
problemáticos da T12 e diagnóstico mecanístico da perda de precisão
quantitativa do preditor \(\rho_1\) quando a referência passa do Modelo D para
o Modelo E**.

A base esperada da `main` é:

```text
8cbde23 feat: audit rho1 with Model E sentinels
```

Ao iniciar, a suíte de referência deve apresentar:

```text
328 passed
```

Esta tarefa deve ser executada integralmente e, depois de auditada localmente,
deve terminar com **commit e push para `origin/main`**, conforme as instruções
da Seção 16.

Não execute T13 nem T14. Em particular, **não avalie nenhum caso novo com
\(N=6\) ou \(N=10\)**.

---

## 1. Motivação e pergunta científica

A T12 auditou 28 sentinelas pré-registrados da calibração \(N\leq4\), usando o
Modelo E como referência completa. A implementação e os diagnósticos numéricos
passaram, mas o gate científico resultou em:

```text
t12_gate_supported=false
NO-GO_T13
```

Três condições falharam:

- convergência confirmada da força de interação em \(22/28=78{,}57\%\);
- RMSE logarítmico da lei congelada igual a aproximadamente \(0{,}781\), acima
  de \(\ln 2\);
- somente \(72{,}73\%\) das previsões aplicáveis dentro de um fator 2.

A T12 também mostrou:

- correlação de Spearman alta, aproximadamente \(0{,}950\);
- nenhum falso seguro no limiar de 5%;
- um falso seguro no limiar de 10%:
  `n2_pair_f0.8_d2.5`;
- contribuição `scattered-scattered` comparável, em vários casos, à correção
  coletiva \(\mathbf F^D-\mathbf F^A\);
- seis casos próximos do contato, todos com \(d_{\min}/a=2.1\), cuja
  convergência da interação não foi confirmada até \(L_{\max}=13\);
- quatro casos adicionais nos quais a interação convergiu, mas o canal
  `scattered-scattered` não.

A T12.1 deve responder, separadamente:

1. **Convergência:** os dez casos problemáticos são apenas casos de convergência
   multipolar lenta ou existe indício de instabilidade/ausência de convergência?
2. **Preditor:** a lei congelada falhou porque o valor de \(\rho_1\) deixou de
   ordenar os erros, porque seus coeficientes foram calibrados contra o alvo
   incompleto do Modelo D, ou porque o novo canal `scattered-scattered`
   introduziu dependência geométrica que \(\rho_1\) sozinho não representa?
3. **Próxima etapa:** os dados justificam uma T12.2 de recalibração controlada
   no domínio \(N\leq4\), ou exigem primeiro um novo descritor físico?

Esta tarefa é **diagnóstica**. Ela não deve forçar a aprovação da T13 e não
deve produzir um novo critério final.

---

## 2. Escopo físico e matemático não negociável

Mantenha exatamente:

- fluido hospedeiro ideal, invíscido, homogêneo, ilimitado e sem perdas;
- acústica linear com dependência temporal \(e^{-i\omega t}\);
- esferas fluidas idênticas, homogêneas, fixas e não sobrepostas;
- onda estacionária nodal;
- centros exatamente em \(z=0\);
- \(a=1\), \(E_0=1\), \(ka=0.1\), \(k=0.1\);
- \(f_0=0\);
- os valores positivos de \(f_1\) e as geometrias congeladas pela T08/T12;
- coeficientes exatos de Mie da T10;
- sistema balanceado de produção da T11.1;
- força multipolar completa e sua decomposição da T11;
- ordenamento modal \(n^2+n+m\);
- redução planar exata \(n+m\) ímpar;
- tolerância de convergência por canal igual a \(10^{-5}\);
- exigência de **duas mudanças sucessivas aplicáveis** abaixo da tolerância.

Não altere:

- os Modelos A, B, C, D ou E;
- a fórmula de Silva–Bruus;
- os coeficientes de Rayleigh ou Mie;
- os operadores de translação;
- a fórmula da força completa;
- a normalização de energia;
- a definição de \(\rho_1\);
- o ajuste congelado da T08;
- os limiares congelados de 1%, 5% e 10%;
- qualquer CSV ou figura das T01–T12.

Não introduza:

- \(ka>0.1\), \(f_0\neq0\) ou contrastes negativos;
- absorção, viscosidade, elasticidade, paredes, streaming ou movimento;
- esferas não idênticas;
- clusters aleatórios ou novas geometrias;
- regularização, clipping, pseudoinversa ou mínimos quadrados no solver;
- extrapolação usada como substituta de uma solução direta;
- ajuste pós-hoc de limiares;
- redes neurais, árvores, polinômios de alta ordem ou busca automática de
  features;
- GMRES, FMM, Numba, GPU ou mudança da arquitetura numérica nesta tarefa.

---

## 3. Auditoria inicial obrigatória

Antes de editar:

1. Leia integralmente:

   - `AGENTS.md`;
   - `README.md`;
   - `TASKS.md`;
   - `docs/CONVENTIONS.md`;
   - `docs/DECISIONS.md`;
   - `docs/HANDOFF.md`;
   - `TAREFA_T08_TRANSFERIBILIDADE_CRITERIO_ACOPLAMENTO.md`;
   - `TAREFA_T09_FUNDAMENTACAO_ANALITICA_RHO1.md`;
   - `TAREFA_T10_COEFICIENTES_EXATOS_MIE.md`;
   - `TAREFA_T11_MODELO_E_REFERENCIA_COMPLETA.md`;
   - `TAREFA_T11_1_ESTABILIZACAO_NUMERICA_MODELO_E.md`;
   - `TAREFA_T12_SENTINELAS_MODELO_E_CRITERIO_RHO1.md`;
   - `src/acoustic_ms/rho_foundation.py`;
   - `src/acoustic_ms/model_e.py`;
   - `src/acoustic_ms/model_e_comparison.py`;
   - `src/acoustic_ms/mie_multiparticle.py`;
   - `src/acoustic_ms/complete_force.py`;
   - `scripts/analyze_t12_model_e_sentinels.py`;
   - os testes das T08, T09, T11, T11.1 e T12.

2. Execute:

   ```bash
   git status --short
   git rev-parse --short HEAD
   git log -1 --oneline
   ```

3. O diretório deve estar limpo e o `HEAD` deve ser `8cbde23`. Se estiver
   apenas desatualizado e sem alterações locais:

   ```bash
   git pull --ff-only origin main
   ```

4. Se houver alterações locais, conflitos ou uma base diferente após o
   `pull`, não descarte nem sobrescreva nada. Pare e informe o estado.

5. Gere um manifesto SHA-256 de todos os arquivos já versionados em:

   ```text
   results/data/
   results/figures/
   ```

   Esse manifesto deve representar os bytes do `HEAD` limpo. A T12.1 pode
   apenas acrescentar artefatos cujo nome comece por `t12_1_`.

6. Instale o projeto e execute a linha de base:

   ```bash
   .venv/bin/python -m pip install -e ".[dev,plot]"
   .venv/bin/python -m pytest -q -W error
   ```

   O resultado esperado é:

   ```text
   328 passed
   ```

7. Confirme que a suíte não reescreveu nenhum artefato versionado. Se a linha
   de base falhar ou algum arquivo antigo mudar, diagnostique e pare antes de
   implementar.

---

## 4. Parte A — Extensão pré-registrada da convergência

### 4.1 Casos autorizados

Não repita a campanha de 28 casos. Estenda **somente** os dez casos abaixo.

#### Seis casos sem convergência confirmada da interação

```text
n2_pair_f1.0_d2.1
n3_compact_f0.8_d2.1
n3_irregular_f1.0_d2.1
n3_linear_f1.0_d2.1
n4_irregular_f0.8_d2.1
n4_linear_f0.8_d2.1
```

#### Quatro casos com interação confirmada, mas `scattered-scattered` não
confirmado

```text
n3_compact_f0.1_d2.1
n4_compact_f0.1_d2.1
n4_irregular_f0.1_d2.1
n4_linear_f0.1_d2.1
```

O conjunto de IDs deve ser validado diretamente contra
`results/data/t12_model_comparison.csv`. Se não coincidir exatamente com o
estado publicado da T12, pare e informe.

### 4.2 Ordens e regra de parada

Para cada um dos dez casos:

1. Reaproveite como proveniência as linhas \(L_{\max}=2,\ldots,13\) já
   existentes em `results/data/t12_model_e_convergence.csv`.
2. Não recalcule nem reescreva essas linhas.
3. Calcule diretamente:

   \[
   L_{\max}=14,15,\ldots,21.
   \]

4. Pare um caso antes de 21 somente depois de confirmar, com duas mudanças
   sucessivas aplicáveis menores ou iguais a \(10^{-5}\), **todos** os quatro
   canais:

   - total;
   - interação;
   - `external-scattered`;
   - `scattered-scattered`.

5. Se algum canal continuar sem confirmação em \(L_{\max}=21\), registre o
   caso como `unconfirmed_at_21`. Não execute \(L_{\max}>21\), não extrapole a
   força e não classifique o caso como divergente.

6. Publique o CSV bruto apenas de forma atômica, depois que todos os dez casos
   tiverem terminado.

### 4.3 Diagnóstico da cauda multipolar

Para cada canal e caso, registre:

\[
\Delta_L
=
\frac{\mathcal R(\mathbf F_L-\mathbf F_{L-1})}
{\mathcal R(\mathbf F_L)},
\]

usando exatamente a regra de escala e aplicabilidade da T12.

Para fins exclusivamente diagnósticos, calcule também:

\[
q_L=\frac{\Delta_L}{\Delta_{L-1}},
\]

quando ambas as quantidades forem positivas, finitas e aplicáveis.

Registre, sem usar para substituir o critério direto:

- mediana de \(q_L\) nas últimas quatro razões aplicáveis;
- mínimo e máximo das mesmas razões;
- presença de oscilação, definida de forma objetiva como pelo menos duas
  alternâncias de crescimento/decrescimento nas últimas cinco mudanças;
- ordem em que cada canal foi confirmado;
- última mudança sucessiva;
- classificação:
  - `directly_confirmed`;
  - `unconfirmed_at_21`;
  - `not_applicable`.

Não ajuste uma lei de convergência e não use Aitken, Richardson ou soma de
série para modificar as forças oficiais.

### 4.4 Diagnósticos numéricos obrigatórios

Em toda nova ordem confirme:

- `production_solver == "balanced_sqrt"`;
- todos os coeficientes e forças finitos;
- \(\kappa_2(A_q)<10\);
- erro retroativo balanceado \(<10^{-12}\);
- fechamentos \(b=a+Ud\) e \(d=Db\) abaixo de \(10^{-12}\);
- resíduo da decomposição da força abaixo de \(10^{-12}\);
- \(\max|F_z|\) compatível com a simetria documentada;
- dimensão do sistema e número de modos ativos coerentes com a base planar.

Convergência lenta com diagnósticos aprovados deve ser descrita como
**truncamento multipolar insuficiente na ordem anterior**, não como
instabilidade do sistema linear.

---

## 5. Parte B — Reconstrução auditada dos 28 sentinelas

Construa uma tabela derivada com uma linha para cada um dos 28 sentinelas:

- para casos não estendidos, use sem recalcular o último estado confirmado da
  T12;
- para os dez casos estendidos, use o último estado diretamente confirmado da
  T12.1;
- se um caso permanecer sem confirmação da interação em \(L_{\max}=21\),
  mantenha suas métricas de erro como não aplicáveis;
- se a interação convergir, mas algum canal de mecanismo não, permita a
  auditoria de erro/limiar e bloqueie apenas a decomposição mecanística.

Recalcule somente os **artefatos derivados da T12.1**. Não modifique
`t12_model_comparison.csv`, `t12_threshold_audit.csv` ou qualquer outro
artefato da T12.

Para cada caso, preserve:

- parâmetros e proveniência da T08/T12;
- ordem final usada;
- flags de convergência por canal;
- forças A, D e E;
- \(\rho_1\), \(\eta\) e \(\Lambda_{\max}\);
- erro da T08 contra D, \(\varepsilon_A^D\);
- erro observado contra E, \(\varepsilon_A^E\);
- previsão congelada da T08;
- razões de aplicabilidade;
- classificações nos limiares congelados.

Reavalie descritivamente o gate da T12 com os casos agora confirmados, mas:

- não altere o registro histórico `NO-GO_T13`;
- não use o resultado para iniciar T13;
- nomeie o novo campo como `t12_1_diagnostic_gate`, não `t12_gate_supported`.

---

## 6. Parte C — Decomposição vetorial do erro

Para cada caso com todos os canais confirmados, defina os campos vetoriais:

\[
\mathbf C_D=\mathbf F^D-\mathbf F^A,
\]

\[
\mathbf C_M=\mathbf F^E_{\mathrm{ext-sc}}-\mathbf F^D,
\]

\[
\mathbf C_S=\mathbf F^E_{\mathrm{ss}},
\]

\[
\mathbf C=\mathbf F^E-\mathbf F^A
=\mathbf C_D+\mathbf C_M+\mathbf C_S.
\]

Use o produto interno global:

\[
\langle\mathbf X,\mathbf Y\rangle
=
\frac1N\sum_{i=1}^{N}
\mathbf X_i\cdot\mathbf Y_i,
\qquad
\mathcal R(\mathbf X)=\sqrt{\langle\mathbf X,\mathbf X\rangle}.
\]

### 6.1 Alinhamentos

Calcule:

\[
\mu_{jk}
=
\frac{\langle\mathbf C_j,\mathbf C_k\rangle}
{\mathcal R(\mathbf C_j)\mathcal R(\mathbf C_k)},
\]

para:

- \((D,M)\);
- \((D,S)\);
- \((M,S)\);
- \((D,C)\);
- \((M,C)\);
- \((S,C)\).

Cada cosseno deve pertencer a \([-1,1]\), salvo tolerância de máquina. Não faça
clipping silencioso; trate apenas desvios compatíveis com arredondamento e
teste-os.

### 6.2 Projeções assinadas

Calcule:

\[
p_j
=
\frac{\langle\mathbf C_j,\mathbf C\rangle}
{\langle\mathbf C,\mathbf C\rangle},
\qquad j\in\{D,M,S\}.
\]

Quando aplicável:

\[
\boxed{p_D+p_M+p_S=1.}
\]

Essas projeções são assinadas e podem ser negativas ou maiores que um. Não as
descreva como percentuais positivos de composição.

### 6.3 Razões de amplitude

Registre:

\[
r_{S/D}=\frac{\mathcal R(\mathbf C_S)}{\mathcal R(\mathbf C_D)},
\qquad
r_{M/D}=\frac{\mathcal R(\mathbf C_M)}{\mathcal R(\mathbf C_D)},
\]

com flags de aplicabilidade e sem piso dimensional artificial.

Mantenha também o `cancellation_ratio` já definido na T12.

### 6.4 Casos que exigem discussão individual

O relatório deve conter uma auditoria específica, baseada nos vetores e não
apenas em amplitudes, de:

```text
n2_pair_f1.0_d6.0
n2_pair_f0.8_d2.5
n2_pair_f1.0_d2.1
```

Explique:

- por que o primeiro produz o maior fator de erro da previsão congelada;
- por que o segundo se torna falso seguro no limiar de 10%;
- se o terceiro preserva o grande erro observado em \(L_{\max}=13\) após a
  extensão direta.

Não atribua causalidade apenas por correlação. Use alinhamentos, projeções e
razões de amplitude para sustentar a interpretação.

---

## 7. Parte D — Diagnóstico controlado do preditor

Esta seção deve distinguir **ordenação**, **calibração** e **informação física
ausente**.

### 7.1 Resíduo da lei congelada

Para cada caso aplicável, defina:

\[
r_{\mathrm{frozen}}
=
\ln\varepsilon_A^E
-\ln\widehat{\varepsilon}_{A,\mathrm{T08}},
\]

com:

\[
\widehat{\varepsilon}_{A,\mathrm{T08}}
=
2.6353684041458636\,
\rho_1^{1.1088518115798773}.
\]

Calcule associações de Spearman entre \(r_{\mathrm{frozen}}\) e:

- \(N\);
- \(f_1\);
- \(d_{\min}/a\);
- \(\rho_1\);
- \(r_{S/D}\);
- \(r_{M/D}\);
- \(\mu_{DS}\);
- \(p_S\);
- `cancellation_ratio`.

As variáveis derivadas do Modelo E são **explicativas**, não preditores
práticos. Rotule-as claramente como `reference_derived`.

### 7.2 Candidatos pré-registrados

Avalie somente os candidatos abaixo. Não acrescente candidatos após observar
os resultados.

#### P0 — lei congelada da T08

Sem qualquer ajuste:

\[
\widehat\varepsilon_{P0}
=
2.6353684041458636\,
\rho_1^{1.1088518115798773}.
\]

#### P1 — potência de \(\eta\)

\[
\widehat\varepsilon_{P1}=C_\eta\,\eta^{p_\eta}.
\]

#### P2 — potência de \(\Lambda_{\max}\)

\[
\widehat\varepsilon_{P2}
=C_\Lambda\,\Lambda_{\max}^{p_\Lambda}.
\]

#### P3 — potência recalibrada de \(\rho_1\)

\[
\widehat\varepsilon_{P3}
=C_\rho\,\rho_1^{p_\rho}.
\]

#### P4 — mudança de alvo a partir do erro contra D

\[
\widehat\varepsilon_{P4}
=Q\,\varepsilon_A^D.
\]

P4 é apenas um diagnóstico: ele requer o Modelo D e não é um substituto barato
para o critério geométrico.

### 7.3 Validação cruzada obrigatória

Para P1–P4, ajuste os coeficientes **somente dentro de cada conjunto de
treino**, em log-espaço, por mínimos quadrados lineares.

Use validação cruzada:

```text
leave-(N,family)-out
```

com os sete estratos da T12:

```text
(2,pair)
(3,compact)
(3,irregular)
(3,linear)
(4,compact)
(4,irregular)
(4,linear)
```

Em cada fold:

- ajuste apenas nos outros seis estratos;
- preveja exclusivamente o estrato removido;
- armazene a previsão out-of-fold de cada caso;
- não use o fold removido para seleção, escala ou tratamento de outlier.

P0 não é reajustado, mas deve ser avaliado sobre o mesmo conjunto aplicável.

Para cada candidato, reporte sobre as previsões out-of-fold:

- número de casos;
- RMSE logarítmico;
- fator mediano;
- percentil 90 do fator;
- fator máximo;
- fração dentro de fator 2;
- Spearman;
- métricas por estrato;
- coeficientes de cada fold.

Também reporte, apenas de forma descritiva, o ajuste global de P1–P4. Não use
o ajuste global como resultado primário.

### 7.4 Proibições contra sobreajuste

Não:

- exclua os três casos especiais;
- aplique pesos escolhidos após ver resíduos;
- crie interceptos específicos por \(N\), família ou caso;
- use interações, splines ou polinômios;
- selecione subconjuntos de dados;
- redefina o erro-alvo;
- altere \(\rho_1\), \(\eta\) ou \(\Lambda_{\max}\);
- derive limiares novos;
- combine features em um novo preditor nesta tarefa.

Se P3 melhorar a calibração, descreva isso como evidência para uma futura
recalibração. Se P3 continuar falhando, descreva isso como evidência de que
\(\rho_1\) isolado pode não carregar informação suficiente.

---

## 8. Critérios de interpretação e recomendação

O relatório deve produzir exatamente uma das recomendações:

### `NEED_MORE_CONVERGENCE`

Use se:

- algum dos 28 casos continuar sem convergência confirmada da interação até
  \(L_{\max}=21\); ou
- qualquer diagnóstico numérico obrigatório falhar.

### `READY_T12_2_RHO1_RECALIBRATION_STUDY`

Use somente se:

- todos os 28 casos tiverem interação confirmada;
- todos os diagnósticos numéricos passarem;
- P3 obtiver, nas previsões out-of-fold:

  \[
  \mathrm{RMSE}_{\log}\leq\ln2
  \]

  e

  \[
  \text{fração dentro de fator 2}\geq0.80;
  \]

- P3 tiver o menor RMSE logarítmico out-of-fold entre P1–P3 ou ficar a no
  máximo \(0.05\) desse menor RMSE;
- a análise mecanística for compatível com uma mudança de calibração, sem
  evidência forte de classes geométricas sistematicamente não representadas.

Essa recomendação autoriza apenas planejar uma T12.2 de calibração mais ampla
em \(N\leq4\). Ela **não autoriza T13**.

### `NEED_NEW_PHYSICS_INFORMED_DESCRIPTOR`

Use se:

- a interação estiver confirmada nos 28 casos e os diagnósticos passarem, mas
  P3 falhar em um dos critérios quantitativos acima; ou
- os resíduos mostrarem dependência sistemática clara do canal
  `scattered-scattered`, de alinhamento/cancelamento ou da família geométrica
  que \(\rho_1\) não representa.

### `INCONCLUSIVE_SMALL_SENTINEL_SET`

Use se os critérios numéricos forem satisfeitos, mas os 28 sentinelas não
permitirem distinguir de forma robusta recalibração de feature ausente.

Não emita `GO_T13` em nenhuma circunstância.

---

## 9. Artefatos obrigatórios

Crie:

```text
src/acoustic_ms/rho1_model_e_diagnostics.py
scripts/analyze_t12_1_rho1_failure.py
tests/test_t12_1_diagnostics.py
tests/test_t12_1_artifacts.py
results/data/t12_1_extended_convergence.csv
results/data/t12_1_convergence_summary.csv
results/data/t12_1_resolved_comparison.csv
results/data/t12_1_mechanism_diagnostics.csv
results/data/t12_1_predictor_diagnostics.csv
results/data/t12_1_out_of_fold_predictions.csv
results/figures/t12_1_rho1_failure_diagnostics.png
TAREFA_T12_1_DIAGNOSTICO_CONVERGENCIA_FALHA_RHO1.md
```

Atualize apenas quando necessário:

```text
src/acoustic_ms/__init__.py
README.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
```

Não crie notebook. Não altere artefatos anteriores.

### 9.1 Conteúdo mínimo dos CSVs

`t12_1_extended_convergence.csv`:

- os dez IDs exatos;
- linhas herdadas \(L=2,\ldots,13\), marcadas como `source=t12`;
- novas linhas, marcadas como `source=t12_1`;
- forças e diagnósticos por ordem;
- mudanças sucessivas, flags e ordens mínimas de confirmação.

`t12_1_convergence_summary.csv`:

- uma linha por caso e canal;
- ordem final;
- ordem de confirmação;
- última mudança;
- estatísticas de \(q_L\);
- flag de oscilação;
- classificação final.

`t12_1_resolved_comparison.csv`:

- exatamente 28 linhas, na ordem do manifesto T12;
- proveniência da força E;
- todas as métricas A–D–E;
- previsões e classificações congeladas;
- aplicabilidade;
- resultado do `t12_1_diagnostic_gate`.

`t12_1_mechanism_diagnostics.csv`:

- \(\mathcal R(\mathbf C_D)\), \(\mathcal R(\mathbf C_M)\),
  \(\mathcal R(\mathbf C_S)\), \(\mathcal R(\mathbf C)\);
- todos os cossenos;
- \(p_D,p_M,p_S\);
- \(r_{S/D}\), \(r_{M/D}\);
- `cancellation_ratio`;
- resíduo do fechamento vetorial;
- flags e razões de não aplicabilidade.

`t12_1_predictor_diagnostics.csv`:

- registros globais out-of-fold;
- registros por estrato;
- coeficientes por fold;
- ajuste global descritivo;
- recomendação final e condições que a sustentam.

`t12_1_out_of_fold_predictions.csv`:

- uma linha por caso e candidato;
- fold;
- alvo;
- previsão;
- resíduo logarítmico;
- fator de erro;
- parâmetros aprendidos sem o fold;
- aplicabilidade.

---

## 10. Figura obrigatória

Gere uma figura determinística 2×2:

1. **Convergência estendida:** \(\Delta_L\) versus \(L_{\max}\) para os seis
   casos sem interação confirmada, com linha horizontal em \(10^{-5}\).
2. **Predição out-of-fold:** \(\varepsilon_A^E\) observado versus previsto
   para P0 e P3, com diagonal e faixas de fator 2.
3. **Mecanismos:** \(r_{S/D}\) versus \(\rho_1\), com cor por família e marcador
   por \(N\).
4. **Resíduo e alinhamento:** \(r_{\mathrm{frozen}}\) versus \(\mu_{DS}\) ou
   \(p_S\), escolhendo entre os dois pela pergunta física mais clara, não pelo
   melhor valor de correlação.

Requisitos:

- escala logarítmica quando matematicamente apropriado;
- rótulos legíveis e sem sobreposição relevante;
- paleta consistente e acessível;
- nenhuma classificação codificada apenas por cor;
- legenda fora dos dados quando necessário;
- margens adequadas;
- nenhum título que reivindique universalidade;
- inspeção visual obrigatória do PNG final.

---

## 11. API e testes

Implemente no módulo importável apenas rotinas científicas reutilizáveis:

- produto interno e projeções de campos vetoriais;
- cossenos com aplicabilidade explícita;
- razões de amplitude com aplicabilidade;
- ajuste log-linear determinístico;
- divisão dos folds `leave-(N,family)-out`;
- métricas out-of-fold;
- diagnóstico da cauda de convergência.

O script deve cuidar de I/O, campanha, agregação e figura.

Teste, no mínimo:

1. validação de shapes, finitude e argumentos;
2. identidade
   \(\mathbf C=\mathbf C_D+\mathbf C_M+\mathbf C_S\);
3. \(p_D+p_M+p_S=1\) quando aplicável;
4. limites dos cossenos;
5. casos ortogonais, alinhados e anti-alinhados construídos analiticamente;
6. razões não aplicáveis sem piso artificial;
7. folds mutuamente exclusivos e cobertura exata dos casos;
8. ausência de vazamento do estrato de teste no ajuste;
9. recuperação de uma potência sintética conhecida;
10. métricas logarítmicas contra um cálculo manual independente;
11. conjunto exato dos dez casos estendidos;
12. conjunto exato dos 28 casos derivados;
13. nenhuma presença de \(N=6\) ou \(N=10\);
14. regra das duas mudanças sucessivas;
15. classificação `unconfirmed_at_21`;
16. preservação dos hashes T01–T12;
17. schema, ordenamento e finitude dos novos CSVs;
18. determinismo dos artefatos derivados;
19. recomendação final coerente com as condições da Seção 8.

Não escreva testes que apenas reproduzam a mesma função de produção sem um
oráculo independente.

---

## 12. Execução reproduzível

O script deve oferecer:

```bash
.venv/bin/python scripts/analyze_t12_1_rho1_failure.py
```

para executar a extensão autorizada e gerar todos os artefatos, e:

```bash
.venv/bin/python scripts/analyze_t12_1_rho1_failure.py --analyze-only
```

para:

- validar a completude do CSV bruto;
- reconstruir os artefatos derivados;
- não executar novas soluções do Modelo E;
- não alterar o CSV bruto.

Depois da campanha:

1. execute `--analyze-only` duas vezes;
2. confirme hashes idênticos de todos os artefatos derivados;
3. confirme que o CSV bruto não mudou;
4. execute a suíte completa com warnings tratados como erros.

Não exija identidade binária entre ambientes diferentes para números de
condição legados mal condicionados. Exija determinismo no mesmo ambiente e
estabilidade física nos campos de força.

---

## 13. Critérios de aceite

A T12.1 só pode ser considerada tecnicamente concluída se:

- a base inicial tiver `328 passed`;
- somente os dez casos autorizados forem estendidos;
- nenhum caso \(N=6,10\) for avaliado;
- nenhuma ordem acima de \(L_{\max}=21\) for executada;
- todos os novos solves usarem `balanced_sqrt`;
- todos os diagnósticos numéricos passarem;
- casos não confirmados forem mantidos como `unconfirmed_at_21`;
- nenhuma força extrapolada for usada;
- os 28 casos forem reconstruídos sem reescrever T12;
- a validação cruzada não tiver vazamento;
- P0–P4 forem os únicos candidatos avaliados;
- os artefatos T01–T12 permanecerem byte-idênticos;
- todos os novos testes passarem com `-W error`;
- `git diff --check` passar;
- os artefatos derivados forem determinísticos no mesmo ambiente;
- a figura for inspecionada visualmente;
- o relatório separar claramente resultado numérico, resultado científico,
  limitações e recomendação;
- a recomendação for exatamente uma das quatro previstas na Seção 8;
- T13 e T14 permanecerem não iniciadas.

O fracasso de um gate científico não significa fracasso de implementação.
Registre o resultado verdadeiro sem modificar critérios para obter aprovação.

---

## 14. Documentação científica obrigatória

Em `TAREFA_T12_1_DIAGNOSTICO_CONVERGENCIA_FALHA_RHO1.md` e
`docs/HANDOFF.md`, documente:

- base Git e ambiente;
- arquivos criados e modificados;
- dez casos estendidos;
- ordens realmente calculadas;
- convergência por canal;
- máximos dos diagnósticos numéricos;
- resultado reconstruído dos limiares congelados;
- métricas P0–P4, com distinção entre in-sample e out-of-fold;
- interpretação das projeções e alinhamentos;
- análise individual dos três dímeros especificados;
- recomendação final;
- limitações:
  - somente 28 sentinelas;
  - \(N\leq4\);
  - \(ka=0.1\);
  - \(f_0=0\);
  - contrastes positivos amostrados;
  - plano nodal;
  - esferas idênticas;
  - ausência de validação externa \(N=6,10\).

Registre explicitamente:

\[
\boxed{
\text{diagnóstico em 28 sentinelas}
\ne
\text{recalibração completa}
\ne
\text{validação externa}
\ne
\text{critério universal}.
}
\]

Atualize `TASKS.md` marcando apenas T12.1 como concluída. Não marque T13,
T14 ou uma eventual T12.2.

---

## 15. Auditoria final antes do Git

Execute:

```bash
.venv/bin/python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
git diff --name-only
```

Além disso:

1. compare o manifesto final dos artefatos T01–T12 com o manifesto inicial;
2. mostre os hashes dos novos artefatos `t12_1_*`;
3. confirme que nenhum arquivo temporário, cache ou ambiente virtual será
   versionado;
4. confirme que não existem resultados de \(N=6\) ou \(N=10\);
5. faça uma auditoria breve do diff;
6. inspecione visualmente a figura;
7. execute duas reconstruções `--analyze-only` e compare os hashes.

Não use:

```bash
git add .
```

Adicione ao staging somente os arquivos pertencentes à T12.1 e as
atualizações documentais autorizadas.

Se este arquivo de prompt estiver localizado dentro do repositório, não o
adicione automaticamente ao commit; o documento oficial da tarefa é
`TAREFA_T12_1_DIAGNOSTICO_CONVERGENCIA_FALHA_RHO1.md`.

---

## 16. Commit e push obrigatórios

Depois que todos os critérios técnicos forem auditados e a conclusão
científica verdadeira estiver documentada, faça o commit:

```bash
git commit -m "feat: diagnose Model E convergence and rho1 failure"
```

Em seguida, envie para a `main` remota:

```bash
git push origin main
```

Se o push for rejeitado por avanço da `main`, não use `--force`. Execute:

```bash
git status --short
git fetch origin
git log --oneline --decorate --graph --max-count=12 --all
```

e pare para informar o conflito.

Ao final, mostre:

```bash
git log -1 --oneline
git status --short
git rev-parse --short HEAD
git rev-parse --short origin/main
```

O estado esperado é:

- commit isolado da T12.1;
- árvore de trabalho limpa;
- `HEAD` e `origin/main` no mesmo commit;
- nenhum uso de force-push;
- nenhum arquivo de `tmp/`, cache, `.venv/` ou artefato estranho incluído.

No relatório final ao usuário, informe de forma concisa:

- hash e mensagem do commit;
- resultado do push;
- total de testes;
- convergência dos dez casos;
- métricas principais de P0–P4;
- recomendação exata da Seção 8;
- lista de arquivos versionados;
- qualquer ressalva numérica ou científica remanescente.
