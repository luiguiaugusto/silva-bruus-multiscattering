# P1.3 — protocolo de validação física de `B_E`

## Estado e separação da campanha

Este documento é o registro prévio e imutável do protocolo da P1.3. A seção
"Protocolo congelado" foi escrita antes da primeira execução do solver Model E
real nesta etapa. A seção "Evidência observada" será preenchida depois, sem
alterar casos, tolerâncias, normalizações ou critérios de aprovação.

Todos os rótulos `DEV-*` abaixo são identificadores locais de
desenvolvimento/oráculo. Eles não são IDs da campanha confirmatória, não foram
derivados dos 102 casos congelados e nenhum valor de resposta dessa campanha é
usado aqui. Nenhum manifesto é habilitado e nenhum artefato é escrito em
`results/` ou `papers/`.

## Protocolo congelado

### Implementações sob teste e oráculos

- Sistema sob teste: `acoustic_ms.model_be.solve_model_be_nodal`, usando o
  `acoustic_ms.model_e.solve_model_e_nodal` real e convergência independente
  para cada par.
- Identidade de dímero e auditoria de ordem comum: chamada direta ao mesmo
  solver Model E, em ordem fixa. Isso audita a composição de `B_E`, não cria
  outro solver.
- Oráculo analítico: `corrected_nodal_pair_forces`, a fórmula corrigida de
  quinta ordem já existente em `acoustic_ms.corrected_pair`; ela não será
  reimplementada.
- Referência de aproximação independente: `nodal_pair_forces`, o Modelo A já
  existente em `acoustic_ms.silva_bruus`.

As implementações históricas B e `B_L` ficam fora do escopo e não serão
modificadas.

### Convenções e normalizações

Adotam-se `a = 1 m` e `E_0 = 1 J/m^3`; portanto `k = ka/a`. Todas as posições
são cartesianas em metros, estão em `z=0`, e os ângulos são medidos a partir de
`+x`. Para um conjunto de vetores `X_i`, define-se

`R(X) = sqrt(mean_i(||X_i||_2^2))`.

Para duas respostas, `Delta_abs(X,Y) = R(X-Y)` e
`Delta_rel(X,Y) = Delta_abs/max(R(X),R(Y))`, quando o denominador é não nulo.
A escala `S` de um gate é `max(R(X),R(Y))` (ou a maior força pertinente ao
invariante). O orçamento de comparação puramente numérica é

`T_num(S) = (5e-10 + 512*epsilon_64) S`.

O termo `512*epsilon_64` cobre arredondamento, transformações e acumulação; o
termo `5e-10` é uma margem conservadora sobre os gates internos de resíduos
(`1e-12`) e os solves lineares bem condicionados esperados. Não é uma
tolerância percentual física.

Todas as chamadas convergentes usam, sem alteração posterior,
`Lmax=2,...,21`, proibição de parada antes de `Lmax=5`, tolerância de mudança
relativa `1e-5` e duas mudanças finais consecutivas aplicáveis em todos os
canais aplicáveis.

### Casos de desenvolvimento/oráculo

| Rótulo | Parâmetros exatos | Finalidade/reuso |
|---|---|---|
| `DEV-GEO` | `ka=0.04`, `f0=0`, `f1=0.35`, `d/a=2.7`, `theta=0.41`; posições `+-1.35*(cos(theta),sin(theta),0)` | Identidade N=2, física do dímero, janela final e ponto alto da comparação analítica. |
| `DEV-GEO-R` | rotação ativa de `DEV-GEO` por `alpha=0.73` | Covariância por rotação. |
| `DEV-GEO-M` | reflexão ativa de `DEV-GEO` por `diag(-1,1,1)` | Covariância por reflexão. |
| `DEV-EQ30-LOW` | como `DEV-GEO`, mas `ka=0.02` | Segundo ponto da auditoria de truncamento da fórmula de quinta ordem. |
| `DEV-NULL` | `ka=0.07`, `f0=f1=0`, `d/a=3.2`, `theta=0.29` | Canais nulos/não aplicáveis. |
| `DEV-WEAK-1..3` | `ka=0.06`, `f0=0`, `d/a=3.3`, `theta=0.27`, com `f1=(0.32,0.08,0.02)` nessa ordem | Tendência contínua `B_E-A -> 0` por contraste dipolar fraco. |
| `DEV-N3` | `ka=0.075`, `f0=0.12`, `f1=0.37`; posições `(-1.55,-0.45,0)`, `(1.25,-0.55,0)`, `(0.20,2.05,0)` | Permutação e ordem comum em triângulo escaleno não sobreposto. |
| `DEV-N3-P` | `DEV-N3` permutado por índices `(2,0,1)` | Covariância por permutação. |

Os testes reutilizam resultados em escopo de módulo. As três transformações de
`DEV-GEO` e a permutação de `DEV-N3` não são novos pontos científicos: são
imagens dos respectivos casos-base.

### Gates congelados

1. **Elegibilidade e diagnósticos.** Todo par aplicável deve ser elegível; seus
   gates numéricos Model E devem passar. Toda ordem final deve estar em
   `[5,21]`. Falha explícita, não convergência ou gate reprovado reprova P1.3.
2. **Identidade N=2.** Em `DEV-GEO`, a força `B_E` deve concordar com a força de
   interação de uma chamada direta ao Model E na ordem final do ledger:
   `Delta_abs <= T_num(S)`.
3. **Covariâncias.** Para rotação, reflexão e permutação, depois de aplicar a
   mesma transformação à resposta de referência,
   `Delta_abs <= T_num(S)`.
4. **Física do dímero.** Somente em `DEV-GEO` — esferas idênticas, fluido sem
   perdas, domínio físico `f0=0`, `-2<=f1<=1`, `d>2a` — planarity, radialidade e
   ação–reação devem ter resíduo absoluto `<= T_num(S)`. Não se generaliza
   ação–reação a partículas diferentes, meios dissipativos ou forças totais de
   campo externo.
5. **Canais nulos.** Em `DEV-NULL`, todos os canais de convergência devem ser
   não aplicáveis, a força deve ser exatamente nula, o caso deve ser elegível e
   a ordem final deve ser 5. Canais nulos/não aplicáveis são dispensados da
   confirmação.
6. **Janela final real.** Em todo canal aplicável de `DEV-GEO`, as duas mudanças
   mais recentes na ordem final devem ser aplicáveis e `<=1e-5`. Confirmações
   históricas anteriores podem permanecer em `confirmation_lmax`, mas não
   substituem essa janela final.
7. **Fórmula corrigida de quinta ordem.** Apenas para `f0=0`, `ka in
   {0.04,0.02}`, `d/a=2.7`, `f1=0.35` e esferas idênticas não sobrepostas:
   `Delta_rel <= (ka)^2` nos dois pontos, o erro baixo deve ser menor que o erro
   alto e a razão `erro_baixo/erro_alto` deve pertencer a `[0.15,0.35]`, janela
   em torno da razão assintótica `1/4` esperada para erro dominante `O(ka^2)`.
   A comparação não valida a fórmula fora desse domínio documentado.
8. **Limite de contraste fraco.** Para `DEV-WEAK-1..3`, tanto
   `R(B_E-A)` quanto `Delta_rel(B_E,A)` devem diminuir estritamente com `|f1|`;
   cada queda absoluta deve exceder `T_num` da escala anterior. Não se impõe
   limiar percentual nem se extrapola um valor no limite.
9. **Ordem comum.** Em `DEV-N3`, cada par é recalculado na maior ordem final
   observada e somado na ordem determinística do ledger. Todos os diagnósticos
   comuns devem passar e
   `Delta_abs(independente,comum) <= 10*1e-5*S + 512*epsilon_64*S`. O fator 10
   é um orçamento de truncamento conservador para três caudas convergidas; não
   é usado nos gates geométricos.

Todos os gates acima foram fixados antes de observar respostas P1.3. Eles não
serão afrouxados depois da execução. Uma reprovação será registrada como
limitação ou `NO_GO`, não tratada ajustando a tolerância.

## Evidência observada

### Proveniência e execução

O protocolo foi congelado no commit `7d4287f`, antes de qualquer chamada real
ao solver nesta etapa. A primeira execução efetiva, em 2026-08-19, usou

```text
python -m pytest -W error -q -s tests/test_model_be_physical.py
```

com o `src/` desta worktree explicitamente à frente do editable install. Não
houve desvio dos parâmetros da tabela congelada. O teste focado terminou em
`8 passed, 1 failed` em `61.57 s`; a única falha foi o gate analítico G7.

### Ordens e diagnósticos Model E

Na tabela, `cond`, `eta_b`, `eta_i` e `eta_s` são, respectivamente, o número de
condição balanceado, o erro backward, o fechamento do campo incidente efetivo
e o fechamento de espalhamento. Para casos com mais de um par registra-se o
pior valor (máximo) entre pares. O resíduo de decomposição da força e
`max|F_z|` foram exatamente zero em todos os pares.

| Caso | Ordens finais por ledger | `cond` máx. | `eta_b` máx. | `eta_i` máx. | `eta_s` máx. | Gate numérico |
|---|---:|---:|---:|---:|---:|---|
| `DEV-GEO` | 12 | 1.0355591301298075 | 1.700e-17 | 1.094e-39 | 8.278e-17 | PASS |
| `DEV-GEO-R` | 12 | 1.0355591301298077 | 2.281e-17 | 1.187e-39 | 8.279e-17 | PASS |
| `DEV-GEO-M` | 12 | 1.0355591301298082 | 7.603e-18 | 1.764e-39 | 8.279e-17 | PASS |
| `DEV-EQ30-LOW` | 12 | 1.0355041594420360 | 1.317e-17 | 9.436e-43 | 8.219e-17 | PASS |
| `DEV-NULL` | 5 | 1.0000000000000000 | 0 | 0 | 0 | PASS |
| `DEV-WEAK-1` | 10 | 1.0135604710482970 | 1.250e-17 | 1.644e-30 | 1.679e-16 | PASS |
| `DEV-WEAK-2` | 10 | 1.0033392946563517 | 2.455e-17 | 1.701e-30 | 8.887e-17 | PASS |
| `DEV-WEAK-3` | 10 | 1.0008317593369032 | 2.113e-17 | 9.315e-30 | 8.189e-17 | PASS |
| `DEV-N3` | 11/11/11 | 1.0317431167040123 | 2.157e-17 | 3.868e-32 | 1.635e-16 | PASS |
| `DEV-N3-P` | 11/11/11 | 1.0317431167040123 | 2.157e-17 | 3.868e-32 | 1.635e-16 | PASS |

No `DEV-GEO`, as primeiras confirmações históricas e as duas mudanças da
janela final foram:

| Canal | `confirmation_lmax` | mudanças finais (`L=11`, `L=12`) | Janela final |
|---|---:|---:|---|
| total | 9 | 6.9082255e-8; 1.2661805e-8 | PASS |
| interação | 9 | 6.9082255e-8; 1.2661805e-8 | PASS |
| externo--espalhado | 8 | 1.0198199e-8; 1.8187132e-9 | PASS |
| espalhado--espalhado | 12 | 3.2512807e-6; 5.9864128e-7 | PASS |

As duas últimas mudanças de todos os canais aplicáveis são aplicáveis e
menores que `1e-5`. Em `DEV-NULL`, os quatro canais são não aplicáveis, não
confirmados e exatamente nulos; a dispensa levou corretamente à parada em
`L=5`.

### Erros, tendências e gates

| Gate | Evidência observada | Resultado |
|---|---|---|
| G1 elegibilidade/diagnósticos | Todos os 14 registros de par foram elegíveis; `cond<1.036`, resíduos internos muito abaixo de `1e-12`, decomposição e `F_z` nulos. | PASS |
| G2 identidade N=2 | `Delta_abs=0`, `Delta_rel=0` contra Model E direto em `L=12`. | PASS |
| G3 rotação | `Delta_abs=9.5771819e-17`, `Delta_rel=2.1214999e-15`. | PASS |
| G3 reflexão | `Delta_abs=3.0046292e-17`, `Delta_rel=6.6557371e-16`. | PASS |
| G3 permutação | `Delta_abs=6.9388939e-18`, `Delta_rel=1.0202757e-16`. | PASS |
| G4 planaridade | resíduo `0`, com escala de força `4.5143448e-2 N`. | PASS |
| G4 radialidade | resíduo `2.3658474e-17 N`; `T_num=2.2577e-11 N`. | PASS |
| G4 ação--reação | resíduo `1.8683544e-17 N`; mesmo `T_num`. | PASS |
| G5 canais nulos | força exatamente nula, quatro canais não aplicáveis, `L=5`. | PASS |
| G6 janela final | valores por canal na tabela anterior; todos simultaneamente confirmados. | PASS |
| G7 fórmula de quinta ordem | ver tabela e análise abaixo. | **FAIL** |
| G8 contraste fraco | correções absoluta e relativa estritamente decrescentes. | PASS |
| G9 ordem comum | ordem comum 11; escala `6.8009991e-2 N`, orçamento `6.8009991e-6 N`, `Delta_abs=Delta_rel=0`; os três diagnósticos passaram. | PASS |

Comparação com `corrected_nodal_pair_forces`, sem reimplementação:

| `ka` | erro absoluto (N) | `Delta_rel` | envelope congelado `(ka)^2` |
|---:|---:|---:|---:|
| 0.04 | 1.1116140914e-3 | 2.4624040583e-2 | 1.6e-3 |
| 0.02 | 1.1169435307e-3 | 2.4774407641e-2 | 4.0e-4 |

O erro baixo não diminuiu: a razão `erro_baixo/erro_alto` é aproximadamente
`1.0061`, fora de `[0.15,0.35]`, e ambos excedem o envelope congelado. O gate
não foi ajustado. A documentação histórica explica que a Eq. (30) corresponde
à redução Rayleigh estrita de ordens ímpares, enquanto Model D planar geral —
e, por extensão, o Model E completo usado aqui — retém canais permitidos
adicionais. Assim, o desvio quase constante não identifica instabilidade de
`B_E`; ele mostra que a hipótese `O((ka)^2)` congelada não é válida para a
comparação entre esses dois espaços de modelo. A fórmula continua sendo um
oráculo truncado útil somente em seu domínio declarado, não uma igualdade para
o Model E completo.

Para o limite de contraste fraco:

| `f1` | `R(B_E-A)` (N) | `Delta_rel(B_E,A)` | escala (N) |
|---:|---:|---:|---:|
| 0.32 | 2.2431083962e-4 | 1.3506513677e-2 | 1.6607604669e-2 |
| 0.08 | 2.4398956585e-6 | 2.3771490039e-3 | 1.0263957600e-3 |
| 0.02 | 2.4433017457e-8 | 3.8178235304e-4 | 6.3997241523e-5 |

Não foi imposto erro percentual final; somente a tendência estrita
pré-registrada foi testada.

### Limitações e decisão

- As afirmações geométricas valem apenas para esferas idênticas, fluido sem
  perdas, plano nodal, domínio de não sobreposição e os casos pequenos
  declarados; não são um teorema para materiais distintos ou forças totais.
- A auditoria de ordem comum foi executada como especificado, mas os três pares
  de `DEV-N3` convergiram todos em `L=11`; por isso o erro comum foi exatamente
  zero e o caso não exercitou ordens finais distintas.
- A tendência de contraste fraco contém três pontos e não constitui ajuste de
  lei de potência nem extrapolação quantitativa ao limite.
- O gate analítico pré-registrado falhou e não pode ser reparado por ajuste
  pós-resposta. Uma nova comparação teria de ser pré-registrada explicitamente
  para o espaço estrito de ordens ímpares, ou então tratar a Eq. (30) apenas
  como diagnóstico com erro de truncamento de canais claramente separado.

Decisão desta P1.3: **`NO_GO_P1.4`** até auditoria científica do gate G7 e da
limitação de ordem comum. Nenhuma atividade da P1.4 foi iniciada.

### Verificação final

- Focado final: `8 passed, 1 failed` em `64.83 s`; falha única em G7.
- Controle histórico isolado após disponibilizar `.venv` na worktree: `1
  passed` em `4.37 s`.
- Suíte completa final com warnings como erros: `534 passed, 1 failed` em
  `141.47 s`; falha única em G7.

A primeira tentativa completa também encontrou a ausência puramente ambiental
de `.venv` na worktree temporária (`533 passed, 2 failed`). Um symlink ignorado
para o ambiente virtual existente removeu essa falha; a repetição isolada e a
suíte final demonstram que ela não era regressão do repositório. Nenhum warning
foi emitido e nenhum arquivo versionado em `results/` ou `papers/` mudou.

## P1.3a — pré-registro da decomposição hierárquica

### Classificação e proveniência

Classificação: **`development_followup`**. Este protocolo é posterior à
observação da falha G7 e, portanto, não é response-blind. Ele foi escrito e
publicado antes de qualquer novo solve da P1.3a. Somente a documentação e
helpers/testes de auditoria podem mudar; APIs de produção permanecerão
inalteradas.

O G7 original permanece congelado sem qualquer alteração retroativa:

| Parâmetro/resultado original | Valor preservado |
|---|---:|
| `d/a`, `f0`, `f1` | `2.7`, `0`, `0.35` |
| `ka` | `0.04`, `0.02` |
| erro relativo observado | `2.4624040583e-2`, `2.4774407641e-2` |
| envelope congelado | `(ka)^2`: `1.6e-3`, `4.0e-4` |
| razão baixo/alto | aproximadamente `1.0061`, fora de `[0.15,0.35]` |
| decisão | **FAIL / `NO_GO_P1.4`** |

Nenhuma tolerância desse gate será alterada. Caso a P1.3a demonstre que ele
compara objetos não equivalentes, o teste será apenas classificado como
`xfail(strict=True)`, mantendo parâmetros e asserts originais visíveis.

### Casos e orientação

Os quatro casos `development_followup` usam `a=1 m`, `E_0=1 J/m^3`,
`d/a=2.7`, `f0=0`, `f1=0.35` e `ka=(0.08,0.04,0.02,0.01)`. Para tornar a
seleção de `m` inequívoca, os centros ficam em `(-1.35,0,0)` e
`(1.35,0,0)`. O escalar reportado é a força na partícula esquerda projetada
na direção da partícula direita; atração é positiva. A covariância rotacional
já auditada em G3 liga esta representação ao dímero oblíquo original, sem
alterar sua força radial.

Nenhum desses rótulos ou resultados pertence aos 102 casos confirmatórios.

### Identidade de reflexão e espaços modais

Será verificada diretamente, com `spherical_harmonic`, a identidade

`Y_l^m(pi-theta,phi)=(-1)^(l+m)Y_l^m(theta,phi)`

para todo `0<=l<=5`, `-l<=m<=l` e os pontos
`(theta,phi)=(0.37,0.23),(0.91,1.17),(1.42,2.03)`. O orçamento absoluto é
`512*epsilon_64*max(1,|Y_l^m|)`.

Em `L=5` ficam congelados:

- espaço planar completo `P={(l,m): l+m ímpar}`;
- espaço reduzido do artigo `R={(l,m): l ímpar, m par}`;
- complemento permitido removido por `R`,
  `Q={(2,+-1),(4,+-1),(4,+-3)}`.

O helper de auditoria montará o mesmo sistema balanceado por translações já
existente, apenas selecionando `P` ou `R`, e reconstruirá BSCs completos com
zeros fora do espaço. Não será criada API de produção. Para `R`, a solução
global de duas partículas também será comparada ao ramo reduzido de uma
partícula já usado por `test_t07_article_formula.py`.

Para cada `ka` serão registrados, por espaço e modelo de coeficiente:

- norma L2 dos coeficientes espalhados totais;
- norma L2 dos coeficientes no setor `Q` e sua razão com a norma total;
- força de interação, externo--espalhado e espalhado--espalhado;
- contribuição modal definida como `F(P)-F(R)` no mesmo observável,
  `L` e coeficiente. Essa diferença inclui corretamente realimentação e termos
  cruzados; não será apresentada como soma independente de forças de modos.

### Pontes congeladas

Cada linha da cadeia altera somente um aspecto. Serão guardados os valores
assinados e, entre cada par, erro absoluto, erro relativo simétrico e ordem
empírica em `ka`:

1. `E_conv:interaction -> E_conv:external_scattered`: isola o canal
   espalhado--espalhado.
2. `E_conv:external_scattered -> E_L5:Mie:P:external_scattered`: isola a ordem
   convergida versus `L=5`.
3. `E_L5:Mie:P -> E_L5:Mie:R`: isola o complemento modal `Q` com coeficientes
   Mie exatos.
4. `E_L5:Mie:R -> E_L5:Rayleigh:R`: isola coeficientes Mie versus Apêndice A
   no mesmo espaço, ordem e funcional de força completa.
5. `E_L5:Rayleigh:R:complete -> E_L5:Rayleigh:R:Appendix-B`: isola a
   aproximação do funcional externo--espalhado usada na redução do artigo.
6. `E_L5:Rayleigh:R:Appendix-B -> Eq.(30)`: isola a expressão fechada e sua
   truncagem em `ka`.

Em paralelo serão calculadas as pontes `Mie:P -> Rayleigh:P` e
`Rayleigh:P -> Rayleigh:R`, tanto para força completa quanto para o funcional
externo--espalhado do Modelo D, evitando atribuir ao coeficiente uma diferença
que venha do espaço modal ou do observável.

A soma assinada das seis diferenças deve reconstruir
`E_conv:interaction-Eq.(30)` dentro de `T_num`. Isso é um gate de fechamento,
não um ajuste dos tamanhos individuais.

### Métricas e critérios congelados

Mantêm-se `R`, `Delta_abs`, `Delta_rel` e
`T_num(S)=(5e-10+512*epsilon_64)S` definidos no protocolo original. Para uma
sequência positiva de erros em halvings de `ka`, a ordem local é

`p_i=log(error(ka_i)/error(ka_i/2))/log(2)`.

Se um erro estiver abaixo de `T_num`, sua ordem será marcada não aplicável, não
forçada a zero. As tolerâncias abaixo são fixadas antes da resposta P1.3a:

1. **Harmônicos.** Todas as identidades de reflexão devem respeitar o orçamento
   de `512 epsilon` acima.
2. **Reprodução do solver.** O helper `Mie:P:L5` deve reproduzir cada canal do
   Model E de produção em `L=5` com `Delta_abs<=T_num`. A identidade
   `interaction=external_scattered+scattered_scattered` deve obedecer o mesmo
   orçamento.
3. **Apêndice A.** Para cada `l=1,...,5`, o erro complexo relativo entre
   coeficientes Mie e Rayleigh deve diminuir e a inclinação de mínimos
   quadrados nos quatro `ka` deve satisfazer `|p-2|<=0.05`, reutilizando o
   critério assintótico já estabelecido em T10.
4. **Redução do artigo.** No espaço `R`, o helper global e o ramo reduzido
   independente do artigo devem concordar em seu funcional comum dentro de
   `T_num`.
5. **Truncagem da Eq. (30).** Somente após espaço, `L`, coeficientes e
   observável estarem equivalentes, o erro do ramo reduzido para a Eq. (30)
   deve diminuir; cada ordem local aplicável deve pertencer a `[1.7,2.3]` e a
   inclinação global deve pertencer ao mesmo intervalo. Nenhum limite
   percentual é imposto.
6. **Paridade física.** A redução para apenas `l` ímpar será considerada
   inválida para o Model E planar completo se o setor `Q` tiver norma resolvida
   (`>512 epsilon` da norma total) e `|F(P)-F(R)|>T_num` nos pontos
   `ka=0.04,0.02`. Se ambos forem nulos dentro desses orçamentos, a redução será
   considerada válida nesses casos. Não se presume o resultado.
7. **Fechamento hierárquico.** A cadeia telescópica deve reconstruir o G7 em
   cada `ka` com erro `<=T_num`; sinal ou normalização não explicados reprovam
   a P1.3a.
8. **Diagnósticos.** Todos os solves Model E e helpers devem ser finitos, com
   resíduos/fechamentos `<=1e-12` e condição balanceada `<10`.

### Nova auditoria de ordem comum

O caso `DEV-N3-MIXED-ORDER` usa exatamente `a=1 m`, `E_0=1 J/m^3`,
`ka=0.1`, `f0=0`, `f1=1` e posições
`(0,0,0)`, `(2.1,0,0)`, `(0,8,0)`. Ele é de desenvolvimento e não pode
alimentar a campanha.

Cada par converge independentemente com o protocolo original (`L=2,...,21`,
parada mínima 5, tolerância `1e-5`, duas mudanças finais). O gate exige:

- todos os pares elegíveis e pelo menos duas ordens finais distintas;
- recálculo dos três pares na maior ordem observada, na ordem do ledger;
- diagnósticos comuns aprovados;
- `Delta_abs(independente,comum) <=
  (10*1e-5+512*epsilon_64)S`, exatamente o orçamento já congelado.

### Regra para classificação do G7 e decisão

O G7 original só receberá `xfail(strict=True)` se: identidade de reflexão,
reprodução do solver, equivalência do ramo reduzido, truncagem da Eq. (30) e
fechamento telescópico passarem; e o setor `Q` demonstrar que os objetos do G7
são não equivalentes. A razão científica deverá constar no marcador e neste
documento. Caso contrário, a falha continuará normal.

- `GO_P1.4`: decomposição completa, todos os gates equivalentes e auditoria de
  ordens distintas aprovados.
- `GO_P1.4_WITH_CONDITIONS`: causa física demonstrada, somente limitação não
  crítica claramente delimitada.
- `NO_GO_P1.4`: qualquer sinal, normalização, paridade ou força permanecer sem
  explicação.

### Evidência P1.3a

`PENDING_FIRST_P1_3A_SOLVE`

### Emenda P1.3a-1 — fallback de ordem comum

Classificação: **`development_followup` pós-resposta**. O primeiro solve do
caso `DEV-N3-MIXED-ORDER` mostrou o resultado global inelegível: o par rígido
próximo `(0,1)`, com `d/a=2.1`, não confirmou convergência até `L=21`. Esse
resultado e a reprovação do candidato inicial ficam preservados; não se soma
uma força parcial e não se aumenta `Lmax`.

Antes de qualquer solve alternativo, fica congelado
`DEV-N3-MIXED-ORDER-FALLBACK`, com `a=1 m`, `E_0=1 J/m^3`, posições
`(0,0,0)`, `(2.7,0,0)`, `(0,8,0)`, `ka=0.04`, `f0=0`, `f1=0.35`. Ele reutiliza
o dímero próximo da própria hierarquia P1.3a e acrescenta somente os dois pares
distantes necessários ao teste de agregação.

Aplicam-se sem alteração: convergência independente em `L=2,...,21`, parada
mínima 5, tolerância `1e-5`, duas mudanças finais, pelo menos duas ordens
finais distintas, recálculo de todos os pares na maior ordem, diagnósticos
aprovados e orçamento `(10*1e-5+512*epsilon_64)S`. O fallback não apaga a falha
do candidato rígido e, mesmo se passar, esta limitação impede `GO_P1.4` sem
condições.
