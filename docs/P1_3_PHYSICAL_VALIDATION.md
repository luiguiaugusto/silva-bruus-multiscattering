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

`PENDING_FIRST_EXECUTION`

Esta seção receberá parâmetros efetivamente executados, ordens finais,
resíduos/diagnósticos, erros absolutos e relativos, resultados dos gates e as
limitações físicas. O preenchimento não modifica o protocolo congelado.
