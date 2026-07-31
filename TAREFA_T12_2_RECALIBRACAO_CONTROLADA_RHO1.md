# Prompt para execução da T12.2 — Recalibração controlada de \(\rho_1\) contra o Modelo E

Você está trabalhando no repositório:

```text
silva-bruus-multiscattering
```

Execute a **T12.2 — estudo confirmatório e pré-registrado de recalibração da
lei quantitativa entre \(\rho_1\) e o erro do Modelo A, usando o Modelo E como
referência**.

A base esperada da `main` é:

```text
59e5ee5 feat: diagnose Model E convergence and rho1 failure
```

Ao concluir e validar a tarefa, faça **um commit isolado e push para
`origin/main`**, conforme a seção final deste prompt.

Não execute T13 nem T14. **Não calcule nenhum caso novo com \(N=6\) ou
\(N=10\)**.

---

## 1. Contexto e pergunta científica

A T12 demonstrou que a calibração original da T08,

\[
\widehat\varepsilon_A^{(0)}(\rho_1)
=2.6353684041458636\,\rho_1^{1.1088518115798773},
\]

não prevê com precisão suficiente o erro do Modelo A quando a referência passa
do Modelo D para o Modelo E. A T12.1 mostrou que:

- a força de interação está convergida nos 28/28 sentinelas;
- a falha é científica, não instabilidade do solver;
- \(\rho_1\) continua sendo informativo e ordena bem os casos;
- o candidato P3, a mesma lei de potência recalibrada contra E, justificou um
  estudo confirmatório;
- os casos `n2_pair_f1.0_d2.1` e `n3_irregular_f1.0_d2.1` ainda têm somente o
  canal interno `scattered-scattered` marcado como `unconfirmed_at_21`, embora
  a força total de interação esteja confirmada.

A T12.2 deve responder, sem tentar modelos sucessivos até obter aprovação:

> Mantendo \(\rho_1\) e a forma funcional de lei de potência já
> pré-especificada, uma nova calibração contra o Modelo E generaliza entre as
> famílias disponíveis em \(N\leq4\), é conservadora nos limiares de 1%, 5% e
> 10% e pode ser congelada antes da validação externa em \(N=6,10\)?

Esta tarefa pode produzir dois resultados legítimos:

```text
GO_T13_WITH_RECALIBRATED_RHO1
```

ou

```text
NO_GO_T13_RHO1_NOT_QUANTITATIVE
```

Não force o primeiro resultado.

---

## 2. Escopo não negociável

Use exclusivamente os **28 sentinelas pré-registrados da T12**, com os valores
de força total confirmados pela T12.1. Não acrescente, substitua ou remova
casos.

Mantenha inalterados:

- toda a física e as convenções das T01–T12.1;
- os Modelos A, B, C, D e E;
- a definição analítica de \(\rho_1\);
- os coeficientes de Mie, operadores de translação e fórmula de força;
- \(a=1\), \(E_0=1\), \(ka=0.1\), \(k=0.1\), \(f_0=0\);
- as geometrias, famílias, contrastes e separações dos sentinelas;
- o erro-alvo \(\varepsilon_A^E\) e o tratamento de escala/logaritmo já
  documentados na T12.1;
- todos os artefatos e resultados das tarefas anteriores.

Não faça nesta tarefa:

- novos solves multipolares, salvo uma regeneração estritamente necessária
  para testar a pipeline existente;
- novos casos, novas famílias ou novos parâmetros físicos;
- seleção automática de variáveis;
- polinômios, splines, redes neurais, árvores, ensembles ou modelos por família;
- clipping de previsões para fazê-las passar;
- exclusão de outliers;
- ajuste separado por \(N\), geometria ou \(f_1\);
- uso de P4 (erro do Modelo D) como preditor final;
- recalibração de \(\rho_1\) em si;
- alteração retroativa da T08, T09, T12 ou T12.1;
- T13 ou T14.

Os dois canais `scattered-scattered` ainda não confirmados devem permanecer
explicitamente sinalizados. Como a força total convergiu, os dois casos entram
na regressão; não declare, porém, que todos os canais internos convergiram.

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
   - as tarefas T08, T09, T11, T11.1, T12 e T12.1;
   - os módulos, scripts, testes e artefatos usados nas análises T08, T12 e
     T12.1.

2. Execute:

   ```bash
   git status --short
   git rev-parse --short HEAD
   git log -1 --oneline
   ```

3. O diretório deve estar limpo e o `HEAD` deve ser `59e5ee5`. Se estiver
   apenas desatualizado e sem mudanças locais, use:

   ```bash
   git pull --ff-only origin main
   ```

4. Se houver alterações locais, conflito, branch incorreta ou base divergente,
   não descarte nada. Pare e informe o estado.

5. Gere um manifesto SHA-256 de todos os artefatos versionados preexistentes em
   `results/data/` e `results/figures/`. A T12.2 só pode adicionar arquivos
   prefixados por `t12_2_`.

6. Instale a base e execute:

   ```bash
   .venv/bin/python -m pip install -e ".[dev,plot]"
   .venv/bin/python -m pytest -q -W error
   ```

   O esperado na base publicada é **354 testes aprovados**. Se a contagem
   diferir apenas por evolução legítima do repositório, registre a explicação;
   se houver qualquer falha, pare antes de implementar.

7. Confirme que os testes não modificaram os artefatos anteriores.

---

## 4. Dados e unidades de validação

### 4.1 Fonte única dos dados

Leia programaticamente os resultados versionados da T12/T12.1; não transcreva
valores manualmente. Para cada sentinela, retenha no mínimo:

- `case_id`;
- \(N\);
- família geométrica;
- \(f_1\);
- \(d_{\min}/a\);
- \(\rho_1\);
- \(\varepsilon_A^E\) da força total;
- status de convergência da força total;
- status dos canais internos.

Faça assertions de que há exatamente 28 `case_id` únicos, somente \(N=2,3,4\),
força total confirmada em 28/28 e os mesmos sete estratos
`(N, family)` usados na T12.

### 4.2 Unidade de generalização

Considere cada combinação `(N, family)` como um grupo. Devem existir exatamente
sete grupos. Toda avaliação principal deve usar **leave-one-group-out** (LOGO):
em cada fold, um grupo inteiro fica fora do ajuste.

É proibido repartir linhas aleatoriamente ou permitir que casos do mesmo grupo
apareçam simultaneamente em treino e teste na métrica confirmatória.

---

## 5. Modelo único pré-especificado

O único modelo candidato confirmatório é:

\[
\boxed{
\log \widehat\varepsilon_A^E
=\beta_0+\beta_1\log\rho_1
}
\]

ou, equivalentemente,

\[
\boxed{
\widehat\varepsilon_A^E=C_E\rho_1^{\alpha_E},
\qquad C_E=e^{\beta_0},\quad\alpha_E=\beta_1.
}
\]

Use regressão linear ordinária não ponderada em espaço logarítmico, reproduzível
com NumPy. Não ajuste interceptos, pesos ou expoentes por grupo.

Use exatamente a mesma definição de \(\varepsilon_A^E\), logaritmo e eventual
`epsilon_floor` adotada e testada na T12.1. Leia esse valor do código ou da
configuração existente; não invente outro. Registre-o nos metadados.

O modelo congelado da T08, P0, é apenas baseline e não pode ser reestimado:

\[
\widehat\varepsilon_A^{(0)}
=2.6353684041458636\rho_1^{1.1088518115798773}.
\]

---

## 6. Protocolo confirmatório LOGO

Para cada um dos sete folds:

1. retire integralmente um grupo `(N, family)`;
2. ajuste \(\beta_0,\beta_1\) somente nos seis grupos restantes;
3. preveja os casos do grupo retirado;
4. salve os coeficientes do fold e as previsões fora da amostra;
5. calcule, também fora da amostra, a previsão P0 nos mesmos casos.

Depois concatene as 28 previsões OOF, preservando a ordem canônica dos
sentinelas.

Calcule para P0 e para o candidato recalibrado:

- RMSE de \(\log\widehat\varepsilon-\log\varepsilon\);
- MAE logarítmico;
- mediana de \(|\log(\widehat\varepsilon/\varepsilon)|\);
- fração dentro de fator 2;
- fração dentro de fator 1,5;
- Spearman entre erro previsto e observado;
- erro máximo de subestimação em escala logarítmica;
- métricas por grupo, apenas como diagnóstico.

Não apresente métricas in-sample como evidência de generalização.

---

## 7. Limiar quantitativo e falsos seguros

Para cada fold e cada tolerância

\[
\tau\in\{0.01,0.05,0.10\},
\]

derive o limiar **somente com os coeficientes de treino daquele fold**:

\[
\rho_{1,\tau}^{(-g)}
=\left(\frac{\tau}{C_E^{(-g)}}\right)^{1/\alpha_E^{(-g)}}.
\]

Esse cálculo só é válido se \(C_E^{(-g)}>0\) e \(\alpha_E^{(-g)}>0\). Caso
contrário, o gate falha.

No grupo retido, classifique como `predicted_safe` quando
\(\rho_1\leq\rho_{1,\tau}^{(-g)}\). Conte:

- número de casos previstos como seguros;
- verdadeiros seguros;
- falsos seguros;
- falsos inseguros;
- pior excesso \(\varepsilon_A^E-\tau\) entre falsos seguros.

Evite um teste vazio: para cada \(\tau\), deve haver ao menos **3 casos OOF
previstos como seguros**, distribuídos em ao menos **2 grupos**. Se não houver,
o critério é inconclusivo para aquela tolerância e o gate final falha.

Não adicione margem conservadora, quantil residual ou fator de segurança depois
de observar os resultados. Uma margem poderá ser estudada futuramente somente
se esta calibração simples falhar e isso for decidido em nova tarefa.

---

## 8. Ajuste final e incerteza descritiva

Somente após produzir e congelar todas as previsões OOF:

1. ajuste a mesma lei nos 28 casos;
2. reporte \(C_E\), \(\alpha_E\) e os limiares finais de 1%, 5% e 10%;
3. reporte a variação dos coeficientes entre os sete folds LOGO;
4. produza intervalos bootstrap de 95% para \(C_E\), \(\alpha_E\) e os três
   limiares, reamostrando **grupos inteiros**, não linhas individuais.

Use semente fixa e pelo menos 10.000 reamostragens válidas. Se uma amostra não
contiver variação suficiente em \(\log\rho_1\) ou produzir coeficiente não
finito, descarte-a e registre quantas tentativas foram necessárias. O bootstrap
é descritivo e não altera o gate.

O ajuste final só pode ser denominado **calibração candidata congelada para
validação externa**. Ele não é uma validação independente, pois usa os 28 casos.

---

## 9. Gate pré-registrado para a T13

Emita `GO_T13_WITH_RECALIBRATED_RHO1` somente se **todos** os itens abaixo forem
satisfeitos pelas previsões OOF do candidato:

1. os 28 casos têm previsão OOF finita e positiva;
2. \(\alpha_E^{(-g)}>0\) e \(C_E^{(-g)}>0\) nos sete folds;
3. RMSE logarítmico \(\leq\ln 2\);
4. pelo menos 85% das previsões ficam dentro de fator 2;
5. Spearman \(\geq0.90\);
6. RMSE logarítmico menor que o P0 OOF;
7. fração dentro de fator 2 maior que a do P0 OOF;
8. **zero falsos seguros** em cada um dos limiares de 1%, 5% e 10%;
9. cobertura mínima não vazia da seção 7 em cada tolerância;
10. nenhum teste científico, numérico ou de integridade falha.

Caso qualquer item falhe, emita:

```text
NO_GO_T13_RHO1_NOT_QUANTITATIVE
```

Nesse caso, conclua explicitamente que \(\rho_1\) pode continuar útil como
indicador ordinal ou mecanístico, mas a lei de potência simples não foi
validada como critério quantitativo autônomo. Não tente outro modelo nesta
tarefa.

O gate não deve ser relaxado depois de conhecidos os resultados.

---

## 10. Implementação e testes

Adapte-se à arquitetura existente e evite duplicar utilitários da T12.1.
Implemente funções puras e testáveis para:

- ajuste da lei log-log;
- geração determinística dos folds LOGO;
- previsões OOF;
- cálculo das métricas;
- inversão da lei para os limiares;
- matriz de classificação segura/insegura;
- bootstrap agrupado;
- avaliação do gate.

Inclua testes, no mínimo, para:

1. recuperação de \(C\) e \(\alpha\) em dados sintéticos exatos;
2. ausência de vazamento entre treino e teste em todos os folds;
3. cada caso aparece exatamente uma vez no teste OOF;
4. invariância à ordem das linhas de entrada;
5. inversão analítica correta dos limiares;
6. rejeição de \(C\leq0\), \(\alpha\leq0\), valores não finitos e duplicatas;
7. contagem correta de falsos seguros e cobertura mínima;
8. métricas conhecidas em exemplo sintético;
9. reprodutibilidade do bootstrap com semente fixa;
10. gate positivo e negativo em fixtures sintéticos;
11. manifesto exato dos 28 sentinelas;
12. preservação dos dois status `scattered-scattered=unconfirmed_at_21`;
13. determinismo dos CSVs e da figura no mesmo ambiente;
14. nenhuma alteração em artefatos anteriores.

Não escreva testes tautológicos que apenas reproduzam o mesmo cálculo da
função sob teste sem um resultado independente ou sintético conhecido.

---

## 11. Artefatos obrigatórios

Crie, com nomes prefixados por `t12_2_`:

### 11.1 Dados

Em `results/data/`:

1. `t12_2_logo_predictions.csv`
   - uma linha por sentinela;
   - grupo retido;
   - valores observado, P0 e recalibrado OOF;
   - razões e resíduos logarítmicos;
   - classificações de segurança nas três tolerâncias;
   - status da força total e dos canais.

2. `t12_2_logo_fits.csv`
   - um registro por fold;
   - grupo retido, tamanhos de treino/teste, \(C_E^{(-g)}\),
     \(\alpha_E^{(-g)}\) e limiares do fold.

3. `t12_2_metrics.csv`
   - métricas globais OOF para P0 e candidato;
   - métricas diagnósticas por grupo claramente rotuladas.

4. `t12_2_safety_audit.csv`
   - tolerância, cobertura, grupos cobertos, verdadeiros/falsos seguros,
     falsos inseguros e pior violação.

5. `t12_2_final_calibration.csv`
   - ajuste nos 28 casos;
   - \(C_E\), \(\alpha_E\), limiares finais;
   - intervalos bootstrap e metadados do bootstrap;
   - valores P0 para comparação.

6. `t12_2_gate.csv`
   - um item por critério;
   - valor observado, limiar, `pass/fail` e justificativa;
   - decisão final inequívoca.

Use ordenação estável, precisão suficiente para round-trip e metadados que
permitam reprodução. Não arredonde antes de calcular métricas ou gate.

### 11.2 Figura

Crie `results/figures/t12_2_rho1_recalibration.png`, em alta resolução, com:

- observado versus previsto OOF, comparando P0 e recalibrado;
- erro observado versus \(\rho_1\), com a curva P0 e a calibração final
  candidata;
- resíduos OOF por grupo;
- auditoria visual dos limiares de 1%, 5% e 10%.

Use escalas logarítmicas quando apropriado, rótulos legíveis e legenda sem
ocultar dados. Marque os dois casos com canal interno não confirmado sem
confundi-los com força total não convergida.

### 11.3 Documentação

Crie `TAREFA_T12_2_RECALIBRACAO_CONTROLADA_RHO1.md` contendo:

- pergunta e protocolo pré-registrado;
- proveniência dos dados;
- equações;
- tabela dos folds e métricas OOF;
- auditoria dos falsos seguros;
- calibração final e incerteza;
- limitações do pequeno conjunto \(N\leq4\);
- decisão do gate;
- consequência exata para a próxima etapa.

Atualize minimamente `README.md`, `TASKS.md`, `docs/DECISIONS.md` e
`docs/HANDOFF.md`. Não reescreva o histórico.

---

## 12. Auditoria final obrigatória

Ao terminar:

1. execute a suíte completa com warnings como erros;
2. execute `git diff --check`;
3. regenere os artefatos T12.2 duas vezes no mesmo ambiente e compare hashes;
4. compare todos os artefatos anteriores com o manifesto inicial;
5. confirme que apenas arquivos T12.2 e as atualizações documentais permitidas
   mudaram;
6. confirme que nenhum caso \(N=6,10\) foi criado ou executado;
7. confira programaticamente todas as condições do gate;
8. inspecione visualmente a figura;
9. registre versões de Python, NumPy e SciPy e explique que igualdade byte a
   byte entre ambientes não é exigida, mas determinismo no mesmo ambiente é.

Se qualquer verificação falhar, corrija dentro do escopo e repita a auditoria.
Não faça commit de uma tarefa com testes falhando.

---

## 13. Commit e push obrigatórios

Somente depois de todas as verificações passarem:

```bash
git status --short
git diff --stat
git diff --check
git add <somente arquivos da T12.2 e atualizações documentais previstas>
git commit -m "feat: recalibrate rho1 against Model E"
git push origin main
```

Regras:

- não use `git add .` sem antes conferir rigorosamente o escopo;
- não inclua arquivos temporários, ambientes virtuais ou caches;
- não altere commits anteriores;
- não use `--amend`, rebase ou `force-push`;
- se o push for rejeitado por avanço remoto, pare e informe; não faça merge ou
  rebase automaticamente;
- se o gate científico falhar, **ainda assim** faça commit e push dos resultados
  válidos com a mesma mensagem: um `NO-GO` é um resultado científico legítimo.

No relatório final ao usuário, informe:

- hash e mensagem do commit;
- confirmação do push;
- contagem de testes;
- arquivos criados/modificados;
- coeficientes e métricas OOF;
- auditoria dos três limiares;
- decisão exata do gate;
- limitações e qualquer ressalva de reprodutibilidade.


---

## Registro da execução

A T12.2 foi executada sobre os 28 sentinelas versionados da T12/T12.1, sem
novos solves e sem acesso ao holdout \(N=6,10\). Os sete folds LOGO retiveram,
em ordem, `n2_pair`, `n3_compact`, `n3_irregular`, `n3_linear`, `n4_compact`,
`n4_irregular` e `n4_linear`; cada fold usou 24 casos de treino e quatro de
teste. Todas as previsões OOF foram finitas e positivas e os coeficientes dos
sete folds foram positivos.

As métricas globais OOF foram:

| modelo | RMSE log | MAE log | dentro de fator 2 | dentro de fator 1,5 | Spearman |
|---|---:|---:|---:|---:|---:|
| P0 congelado | 0.8106033557995027 | 0.6330031366128833 | 0.6071428571428571 | 0.4285714285714286 | 0.9709906951286261 |
| \(\rho_1\) recalibrado | 0.6458489737104012 | 0.3786559305780956 | 0.9285714285714286 | 0.7142857142857143 | 0.9600437876299945 |

A auditoria de segurança OOF encontrou 7, 14 e 20 casos previstos como seguros
nos limiares de 1%, 5% e 10%, cobrindo os sete grupos em cada tolerância. Os
limiares de 1% e 5% tiveram zero falsos seguros. No limiar de 10%,
`n2_pair_f0.8_d2.5` permaneceu falso seguro, com excesso
0.020573189849995427. Esse único resultado faz falhar o critério 8 do gate.

O ajuste descritivo nos 28 casos é

\[
\widehat\varepsilon_A^E
=14.73950709797405\,\rho_1^{1.4226504975598322},
\]

com limiares candidatos 0.005926947606709601, 0.01837157635582504 e
0.029905042165737895 para 1%, 5% e 10%. O bootstrap agrupado usou semente 1202
e 10.000 amostras válidas em 10.000 tentativas. Os intervalos percentis de 95%
foram [7.57978280636931, 47.07370137527077] para o prefator e
[1.24964935549445, 1.728753959643017] para o expoente. O bootstrap é apenas
descritivo e não altera o gate.

A decisão exata é:

```text
NO_GO_T13_RHO1_NOT_QUANTITATIVE
```

Assim, \(\rho_1\) permanece informativo como indicador ordinal e mecanístico,
mas a lei de potência simples não foi validada como critério quantitativo
autônomo. Os dois casos com `scattered-scattered=unconfirmed_at_21` permanecem
sinalizados, embora suas forças de interação estejam diretamente confirmadas.
T13 e T14 não foram iniciadas.
