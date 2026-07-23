# T04 - Forca de interacao por multiplo espalhamento Rayleigh (`Lmax=1`)

## 1. Objetivo

Implementar, testar e documentar o calculo da forca acustica de interacao no
plano nodal a partir dos coeficientes multipolares acoplados produzidos pela
T03.

Esta tarefa deve entregar o primeiro calculo de forca do **Modelo C**:

\[
\boxed{
\text{multiplo espalhamento global com }L_{\max}=1
\text{ e todas as ordens de reespalhamento}
}
\]

O escopo cientifico desta tarefa termina nos benchmarks de **uma e duas
particulas**. A implementacao pode ser naturalmente vetorizada para receber os
coeficientes de \(N\) particulas, mas a T04 nao deve produzir, interpretar ou
publicar resultados de forca para \(N\geq3\). Triades e a contribuicao
irredutivel de tres corpos pertencem a uma tarefa posterior.

Ao final, o repositorio deve ser capaz de:

1. resolver o sistema acoplado Rayleigh da T03;
2. reexpandir, em cada centro, o campo espalhado pelas demais particulas ate
   os modos locais necessarios para a forca;
3. calcular a forca de interacao bidimensional em newtons;
4. recuperar uma referencia analitica independente para duas particulas;
5. recuperar Silva-Bruus nos limites de contraste fraco e acoplamento fraco;
6. manter integralmente aprovadas as T01, T02, T03 e T03.1.

---

## 2. Estado inicial e preparacao obrigatoria

A base cientifica aprovada esta na `main`, no commit:

```text
58039bb752f644753af90e33585a3213519dd59a
```

Antes de alterar qualquer arquivo:

1. execute `git status --short`;
2. execute `git branch --show-current`;
3. execute `git rev-parse HEAD`;
4. preserve qualquer alteracao preexistente do usuario;
5. leia integralmente:
   - `AGENTS.md`;
   - `docs/CONVENTIONS.md`;
   - `docs/DECISIONS.md`;
   - `docs/HANDOFF.md`;
   - a especificacao da T03;
   - os modulos e testes implementados na T03;
6. execute:

   ```bash
   python -m pytest -q -W error
   ```

7. confirme a linha de base de **58 testes aprovados**.

Se a arvore estiver suja, a branch nao for `main`, o `HEAD` nao corresponder ao
estado esperado ou houver divergencia remota que exija sobrescrever trabalho,
pare e relate o problema. Nao use `reset --hard`, `checkout --` ou outro
comando destrutivo.

---

## 3. Regime fisico e convencoes que devem ser preservados

Use exatamente as convencoes ja aprovadas:

- dependencia temporal:

  \[
  p(\mathbf r,t)=\operatorname{Re}
  \left\{p(\mathbf r)e^{-i\omega t}\right\};
  \]

- fluido ideal, inviscido, homogeneo e ilimitado;
- acustica linear e harmonica;
- esferas identicas, compressiveis e nao sobrepostas;
- regime de Rayleigh:

  \[
  0<ka\leq0.1;
  \]

- onda estacionaria normalizada:

  \[
  \widetilde p^{\mathrm{ext}}=\sin(kz);
  \]

- centros no plano nodal \(z=0\);
- energia acustica explicita:

  \[
  E_0=\frac{\rho_0|v_0|^2}{4};
  \]

- vetor entre fonte e prova:

  \[
  \mathbf d_{ij}=\mathbf r_i-\mathbf r_j,
  \qquad
  \widehat{\mathbf d}_{ij}
  =\frac{\mathbf d_{ij}}{d_{ij}},
  \]

  apontando da fonte \(j\) para a particula de prova \(i\);
- forca radial negativa significa atracao;
- nao sobreposicao:

  \[
  d_{ij}\geq2a;
  \]

- base de harmonicos esfericos complexos ortonormais com fase de
  Condon-Shortley;
- ordenacao:

  \[
  \operatorname{index}(\ell,m)=\ell^2+\ell+m;
  \]

- operador de translacao orientado como
  `target <- source`, com linhas correspondentes aos modos do alvo e colunas
  aos modos da fonte.

Nao trate \(ka\), \(kd\) e \(a/d\) como variaveis independentes:

\[
kd=(ka)\frac da.
\]

---

## 4. Distincao obrigatoria: forca total, primaria e de interacao

O artigo de 2026 apresenta a forca total de radiacao na Eq. (21) e a forca de
interacao, na ordem assintotica adotada, na Eq. (22).

A T04 deve implementar a **forca de interacao especializada ao plano nodal**,
isto e, a Eq. (27), e nao uma nova avaliacao irrestrita da Eq. (21).

No plano nodal:

- a forca primaria da particula isolada e nula;
- a T04 deve retornar zero para \(N=1\);
- os coeficientes acoplados da T03 ja incluem todas as cadeias de
  reespalhamento permitidas por \(L_{\max}=1\);
- a expressao da forca mantem os termos cruzados
  campo externo-campo espalhado da Eq. (22);
- nao devem ser acrescentados termos quadraticos
  espalhado-espalhado que nao pertencem a Eq. (22);
- nao se deve chamar o resultado de forca total fora do plano nodal.

Em particular, **nao** substitua simplesmente o campo incidente efetivo
completo na Eq. (21) e mantenha todos os produtos quadraticos. Isso produziria
outro nivel de aproximacao e impediria a comparacao controlada com as
Eqs. (27)-(30) do benchmark de 2026.

---

## 5. Por que a avaliacao da forca exige modos locais `ell=2`

O solver da T03 deve continuar com:

\[
L_{\max}^{\mathrm{scatter}}=1,
\]

isto e, quatro incognitas por particula:

\[
(0,0),\quad(1,-1),\quad(1,0),\quad(1,1).
\]

Entretanto, a forca depende do gradiente da velocidade espalhada, ou,
equivalentemente, de derivadas segundas da pressao. Por isso, o campo
espalhado pelas outras particulas deve ser **avaliado localmente** no centro
da particula-alvo ate \(\ell=2\):

\[
b_{nm}^{(i,\mathrm{sc})}
=
\sum_{j\ne i}
\sum_{n'=0}^{1}\sum_{m'=-n'}^{n'}
S_{n'm',nm}^{(j\rightarrow i)}
s_{n'm'}^{(j)},
\qquad n\leq2.
\]

Use a API retangular que ja existe:

```python
translation_matrix(
    k,
    target_position,
    source_position,
    target_lmax=2,
    source_lmax=1,
)
```

Portanto:

\[
\boxed{
L_{\max}^{\mathrm{scatter}}=1
\quad\text{e}\quad
L_{\max}^{\mathrm{evaluation}}=2
}
\]

nao significam que o solver foi promovido a \(L_{\max}=2\). Os modos locais
\(b_{2m}\) sao apenas observaveis do campo regular reexpandido; nao sao novas
incognitas de espalhamento e nao devem realimentar o sistema da T03.

O campo proprio da particula-alvo deve ser excluido dessa soma. Nunca avalie a
onda de Hankel singular da propria particula em seu centro.

---

## 6. Formula de forca a implementar

Para a configuracao nodal planar da T03, a Eq. (27) do artigo de 2026 e:

\[
\mathbf F_i^{\mathrm{int}}
=
4\pi a^3E_0
\operatorname{Re}
\left[
i f_1^*
(\mathbf e_z\cdot\nabla_i)
\widetilde{\mathbf v}_i^{\mathrm{sc}}
\right]_{\mathbf r_i=0}.
\]

Com

\[
\widetilde{\mathbf v}^{\mathrm{sc}}
=-\frac{i}{k}\nabla\widetilde p^{\mathrm{sc}},
\]

e usando as expansoes cartesianas na origem,

\[
j_2(kr)Y_2^1
=
-\frac{k^2}{15}
\sqrt{\frac{15}{8\pi}}\,
z(x+iy)
+O(r^4),
\]

\[
j_2(kr)Y_2^{-1}
=
\frac{k^2}{15}
\sqrt{\frac{15}{8\pi}}\,
z(x-iy)
+O(r^4),
\]

obtem-se:

\[
F_{i,x}
=
C_F
\operatorname{Re}
\left[
f_1^*
\left(
b_{2,-1}^{(i,\mathrm{sc})}
-b_{2,1}^{(i,\mathrm{sc})}
\right)
\right],
\]

\[
F_{i,y}
=
C_F
\operatorname{Re}
\left[
-i f_1^*
\left(
b_{2,1}^{(i,\mathrm{sc})}
+b_{2,-1}^{(i,\mathrm{sc})}
\right)
\right],
\]

com

\[
\boxed{
C_F=
\frac{\sqrt{30\pi}}{15}
k\,a^3E_0
}.
\]

No escopo planar:

\[
F_{i,z}=0.
\]

No codigo atual, \(f_1\) e real, mas preserve a conjugacao indicada pela
derivacao ao organizar a formula.

Nao use diferencas finitas, quadratura angular ou avaliacao em pontos
deslocados no codigo de producao. A formula deve ser avaliada diretamente
pelos coeficientes \(b_{2,\pm1}\).

---

## 7. Referencia independente para os testes

Os testes nao podem validar a formula apenas repetindo a mesma combinacao de
coeficientes usada no modulo de producao.

No plano nodal e com \(L_{\max}^{\mathrm{scatter}}=1\), apenas o coeficiente
\(s_{10}^{(j)}\) sobrevive. O campo de cada fonte pode ser escrito diretamente
em coordenadas cartesianas como:

\[
\widetilde p_j^{\mathrm{sc}}(\mathbf r)
=
s_{10}^{(j)}
h_1^{(1)}(kR_j)
\sqrt{\frac{3}{4\pi}}
\frac{z}{R_j},
\]

onde:

\[
R_j=|\mathbf r-\mathbf r_j|.
\]

Defina:

\[
Q(x)
=
\frac{h_1^{(1)\prime}(x)}{x}
-\frac{h_1^{(1)}(x)}{x^2}.
\]

Uma expressao cartesiana independente para a forca e:

\[
\boxed{
\mathbf F_i^{\mathrm{int},L=1}
=
4\pi k a^3E_0
\sqrt{\frac{3}{4\pi}}
\operatorname{Re}
\left[
f_1^*
\sum_{j\ne i}
s_{10}^{(j)}
Q(kd_{ij})
\widehat{\mathbf d}_{ij}
\right]
}.
\]

Implemente essa expressao **somente nos testes**, como oraculo independente.
Ela nao deve chamar:

- `translation_matrix`;
- `separation_coefficient`;
- coeficientes de Gaunt;
- a funcao de forca de producao.

Ela pode usar diretamente:

```python
spherical_hankel1(1, x)
spherical_hankel1(1, x, derivative=True)
```

Essa comparacao deve ser feita para pares alinhados e rotacionados, de modo a
testar simultaneamente:

- o sinal;
- as fases de \(b_{2,\pm1}\);
- a reconstrucao das componentes \(x\) e \(y\);
- a orientacao `target <- source`.

---

## 8. Benchmark analitico escalar de duas particulas

Para duas particulas identicas em:

\[
\mathbf r_1=(-d/2,0,0),
\qquad
\mathbf r_2=(d/2,0,0),
\]

o sistema \(L_{\max}=1\) reduz-se a:

\[
s_{10}
=
\frac{s_1a_{10}}
{1-s_1\left[
h_0^{(1)}(kd)+h_2^{(1)}(kd)
\right]},
\]

com:

\[
s_1=\frac{i}{6}f_1(ka)^3,
\qquad
a_{10}=\sqrt{12\pi}.
\]

A forca radial assinada ao longo de
\(\widehat{\mathbf d}_{12}\) deve ser:

\[
\boxed{
F_{\parallel}^{L=1}
=
4\pi k a^3E_0
\sqrt{\frac{3}{4\pi}}
\operatorname{Re}
\left[
f_1^*s_{10}Q(kd)
\right]
}.
\]

Ela deve ser negativa no regime atrativo de curto alcance.

Esse benchmark deve ser calculado no teste a partir das equacoes escalares
acima, sem montar a matriz global \(8\times8\) e sem usar a funcao de forca de
producao.

### Valores de regressao no contato

Para:

\[
a=1,\qquad
E_0=1,\qquad
ka=0.1,\qquad
d/a=2,\qquad
kd=0.2,
\]

os valores radiais assinados do modelo MS Rayleigh \(L_{\max}=1\) sao:

| \(f_1\) | \(F_{\parallel}^{L=1}\) |
|---:|---:|
| 0.1 | \(-0.011936371917121\) |
| 0.4 | \(-0.194729303800953\) |
| 0.8 | \(-0.799842697325624\) |
| 1.0 | \(-1.26676999261163\) |

Como \(a=1\,\mathrm m\) e \(E_0=1\,\mathrm{J\,m^{-3}}\) nessa regressao
normalizada, os valores numericos acima correspondem a newtons. Eles devem ser
usados como regressao de sinal e normalizacao, nao como parametros
experimentais.

Use tolerancia relativa maxima de \(2\times10^{-12}\) para esses valores.

---

## 9. Benchmark direto da Eq. (28) do artigo de 2026

Para o par alinhado ao eixo \(x\), calcule:

\[
b_{21}^{(2\rightarrow1)}
\]

pela reexpansao retangular `target_lmax=2, source_lmax=1`.

Verifique diretamente:

\[
\boxed{
\mathbf F_1^{\mathrm{int}}
=
-\frac{2\sqrt{30\pi}}{15k^2}
(ka)^3E_0
\operatorname{Re}
\left[
f_1^*b_{21}^{(2\rightarrow1)}
\right]\mathbf e_x
}
\]

e:

\[
\mathbf F_2^{\mathrm{int}}=-\mathbf F_1^{\mathrm{int}}.
\]

No par disposto como acima, teste tambem as relacoes de simetria dos
coeficientes locais:

\[
b_{2,-1}^{(1,\mathrm{sc})}
=-b_{2,1}^{(1,\mathrm{sc})},
\]

e que os modos locais proibidos pela simetria nao contribuem acima da
tolerancia numerica.

Esse teste conecta explicitamente a implementacao a Eq. (28), mas nao
substitui o oraculo cartesiano independente da Secao 7.

---

## 10. Relacao esperada com os modelos A e B

A T04 implementa:

\[
\mathbf F^{\mathrm{MS},L=1},
\]

que nao e identica, em geral, a:

\[
\mathbf F^{\mathrm{SB}}
\]

nem a:

\[
\mathbf F^{\mathrm{corr},L=5}.
\]

As comparacoes devem ser interpretadas como:

\[
\mathbf F^{\mathrm{MS},L=1}
-\mathbf F^{\mathrm{SB}}
\quad\Longrightarrow\quad
\text{reespalhamento dentro da base Rayleigh},
\]

enquanto:

\[
\mathbf F^{\mathrm{corr},L=5}
-\mathbf F^{\mathrm{MS},L=1}
\]

tambem contem a correcao multipolar do benchmark de duas particulas.

Portanto:

- nao force igualdade entre a T04 e a formula corrigida da T02;
- nao altere a T02 para faze-la coincidir com \(L_{\max}=1\);
- use a T02 apenas como valor comparativo no CSV;
- exija a recuperacao de Silva-Bruus quando o reespalhamento se torna
  perturbativamente fraco.

Testes obrigatorios desses limites:

1. \(f_1\rightarrow0\), com \(ka\), \(kd\) e \(d/a\) consistentes:

   \[
   \frac{
   |F^{\mathrm{MS},L=1}-F^{\mathrm{SB}}|
   }{
   |F^{\mathrm{MS},L=1}|
   }
   \rightarrow0;
   \]

2. \(a/d\rightarrow0\), longe de um zero da forca:

   \[
   F^{\mathrm{MS},L=1}
   \rightarrow
   F^{\mathrm{SB}}.
   \]

Para o teste de contraste fraco, pode-se usar \(f_1=10^{-7}\) e tolerancia
relativa de \(10^{-8}\). Para o limite de separacao grande, escolha um ponto
longe de zeros e justifique numericamente a tolerancia.

---

## 11. API e arquivos de producao

Crie:

```text
src/acoustic_ms/force.py
```

A API publica deve conter, no minimo, uma funcao de alto nivel com nome
semantico equivalente a:

```python
solve_rayleigh_nodal_interaction_forces(
    positions_xyz,
    k,
    radius,
    energy_density,
    f0,
    f1,
    lmax=1,
)
```

Ela deve:

1. chamar `solve_rayleigh_nodal` exatamente uma vez;
2. construir, para cada alvo, a soma dos campos espalhados por todas as outras
   particulas com:

   ```python
   target_lmax=2
   source_lmax=1
   ```

3. excluir o campo proprio;
4. calcular \(F_x\) e \(F_y\) analiticamente pelos coeficientes
   \(b_{2,\pm1}\);
5. retornar os coeficientes do campo e os observaveis de forca de maneira
   auditavel.

Use uma `dataclass(frozen=True)` com nome semantico equivalente a:

```python
RayleighNodalInteractionResult
```

e campos publicos equivalentes a:

```python
solution: RayleighNodalSolution
local_scattered_coefficients: np.ndarray  # shape (N, 9), ell <= 2
forces_xy: np.ndarray                     # shape (N, 2), SI: N
```

O nome exato pode ser ajustado apenas se houver uma razao clara de consistencia
com a API existente. Nao esconda os coeficientes locais \(b_{nm}\), pois eles
sao necessarios para auditar a Eq. (28) e serao reutilizados quando o
truncamento multipolar for ampliado.

Reutilize:

- `solve_rayleigh_nodal`;
- `translation_matrix`;
- `mode_index`;
- as validacoes ja existentes sempre que adequado.

Evite duplicar o solver, a indexacao ou o operador de translacao em
`force.py`.

Atualize:

```text
src/acoustic_ms/__init__.py
```

para exportar a nova API, sem remover ou renomear simbolos publicos anteriores.

---

## 12. Validacao de entradas e unidades

A nova API deve rejeitar:

- `energy_density < 0`;
- valores nao finitos;
- `k <= 0`;
- `radius <= 0`;
- `ka > 0.1`;
- `f1` fora de \([-2,1]\);
- centros fora do plano nodal;
- centros coincidentes;
- sobreposicao \(d<2a\);
- `lmax != 1`;
- formas de arrays inconsistentes.

As validacoes geometricas e do regime ja fornecidas por
`solve_rayleigh_nodal` devem ser reutilizadas, sem manter duas implementacoes
divergentes.

Confirme dimensionalmente:

\[
[k a^3E_0]
=
\mathrm{m^{-1}m^3J\,m^{-3}}
=
\mathrm N.
\]

Teste:

- linearidade em \(E_0\);
- forca nula para \(E_0=0\);
- forca nula para \(f_1=0\);
- independencia de \(f_0\) na configuracao nodal planar `Lmax=1`;
- escala \(a^2E_0\) quando \(ka\) e \(d/a\) sao mantidos fixos.

---

## 13. Testes obrigatorios

Crie:

```text
tests/test_force.py
```

Inclua, no minimo, os seguintes grupos.

### A. Particula isolada

- `local_scattered_coefficients` deve ser zero;
- `forces_xy` deve ser zero;
- a solucao da T03 deve permanecer identica ao benchmark de particula
  isolada;
- documente que isso corresponde a forca de interacao e a forca primaria
  nulas no centro exato do no.

### B. Oraculo cartesiano independente

Compare a forca de producao com a expressao da Secao 7 para:

- um par alinhado ao eixo \(x\);
- um par em orientacao obliqua;
- pelo menos tres combinacoes de \(ka\), \(d/a\) e \(f_1\);
- um caso proximo ao contato e um caso mais afastado.

Use tolerancia relativa maxima de \(2\times10^{-12}\), salvo justificativa
numerica explicita.

### C. Solucao escalar exata do par

Monte a referencia da Secao 8 diretamente no teste e compare:

- os dois coeficientes \(s_{10}\);
- a forca radial assinada;
- a direcao vetorial;
- os quatro valores de regressao no contato.

Nao use o solver matricial para construir o valor esperado.

### D. Eq. (28)

Verifique diretamente a Eq. (28) com
\(b_{21}^{(2\rightarrow1)}\), incluindo sinal, prefator, direcao e unidades.

### E. Silva-Bruus

Verifique:

- limite \(f_1\rightarrow0\);
- limite \(a/d\rightarrow0\);
- atracao no limite \(kd\ll1\);
- coerencia do sinal radial com a T01.

### F. Simetrias do par

Verifique:

- acao-reacao para o par identico;
- invariancia por translacao comum;
- covariancia por rotacao no plano;
- independencia da ordem das particulas apos a permutacao correspondente dos
  resultados;
- ausencia de componente perpendicular a linha de centros no par.

Nao generalize automaticamente o teste de acao-reacao para clusters
multibody nesta tarefa. Essa propriedade deve ser analisada separadamente
quando os termos de tres corpos forem estudados.

### G. Escalas e dominio

Verifique:

- linearidade em `energy_density`;
- escala \(a^2\) a \(ka\) e \(d/a\) fixos;
- zeros para `energy_density=0` e `f1=0`;
- independencia de `f0`;
- rejeicao de entradas invalidas.

### H. Nao regressao

Todos os 58 testes anteriores devem continuar passando sem alteracao de
tolerancias para esconder falhas.

Nao use:

- `skip`;
- `xfail`;
- supressao de warnings;
- valores esperados gerados pela mesma funcao sob teste.

---

## 14. Script e CSV deterministico

Crie:

```text
scripts/validate_t04_force.py
```

O script deve usar somente a API publica do pacote para gerar:

```text
results/data/t04_pair_force_validation.csv
```

Inclua:

1. uma particula isolada;
2. pares com \(ka=0.1\), \(d/a=2\) e:

   \[
   f_1=0.1,\ 0.4,\ 0.8,\ 1.0;
   \]

3. pelo menos dois pares adicionais com outras separacoes;
4. para cada par:
   - `ka`;
   - `d_over_a`;
   - `kd`;
   - `f0`;
   - `f1`;
   - componentes das forcas nas duas particulas;
   - forca radial MS `Lmax=1`;
   - referencia escalar independente;
   - erro relativo entre as duas;
   - forca Silva-Bruus;
   - forca corrigida T02 `Lmax=5`;
   - residuo do solver;
   - numero de condicao.

Use terminadores de linha `\n` e ordem fixa de colunas e linhas para tornar o
arquivo deterministico.

O script deve imprimir um resumo curto e falhar de forma explicita se o erro
contra a referencia escalar ultrapassar a tolerancia.

Nao gere ainda:

- mapas de erro;
- varreduras grandes;
- figuras de triades;
- resultados para \(N\geq3\);
- trajetorias.

Uma figura nao e exigida na T04.

---

## 15. Documentacao

Atualize somente o necessario.

### `docs/CONVENTIONS.md`

Registre:

- a distincao entre a Eq. (21) e a forca de interacao da Eq. (22);
- a especializacao nodal da Eq. (27);
- \(L_{\max}^{\mathrm{scatter}}=1\) versus
  \(L_{\max}^{\mathrm{evaluation}}=2\);
- a combinacao de \(b_{2,\pm1}\) usada em \(F_x\) e \(F_y\);
- que a saida e uma forca bidimensional em newtons;
- que o campo proprio e excluido.

### `docs/DECISIONS.md`

Registre que a T04:

- implementa o Modelo C no nivel Rayleigh;
- usa a Eq. (22)/(27), sem termos quadraticos espalhado-espalhado;
- nao amplia o solver de producao alem de `Lmax=1`;
- ainda nao apresenta resultados de tres corpos.

### `docs/HANDOFF.md`

Adicione uma secao T04 contendo:

- arquivos alterados;
- equacoes implementadas;
- valores numericos dos benchmarks;
- erros maximos medidos;
- comandos executados;
- contagem final de testes;
- limitacoes restantes.

### `README.md`

Corrija o texto que ainda afirma que multiple scattering esta fora do escopo e
documente brevemente a nova API de forca.

### `src/acoustic_ms/__init__.py`

Atualize o docstring do pacote e exporte a API da T04.

Nao reescreva o historico das tarefas anteriores.

---

## 16. Restricoes

- Nao altere as formulas aprovadas das T01 e T02.
- Nao altere o nucleo cientifico aprovado da T03, salvo se um teste novo
  demonstrar um defeito real e reproduzivel.
- Nao aumente o solver para `Lmax=2`, `3` ou `5`.
- Nao adicione forca primaria fora do plano nodal.
- Nao adicione termos espalhado-espalhado.
- Nao calcule torque.
- Nao calcule trajetorias ou dinamica.
- Nao implemente viscosidade, streaming, paredes ou hidrodinamica.
- Nao produza resultados cientificos para \(N\geq3\).
- Nao implemente a contribuicao irredutivel de tres corpos.
- Nao adicione novas dependencias.
- Nao otimize prematuramente e nao introduza solver esparso, GMRES, Numba ou
  paralelizacao.
- Nao use notebooks como unica implementacao.
- Nao esconda falhas alterando tolerancias sem justificativa.

Se um novo teste revelar um erro real em uma convencao da T03, pare, documente
o contraexemplo minimo e explique o impacto antes de modificar o nucleo
aprovado.

---

## 17. Criterios de aceitacao

A T04 somente esta concluida quando:

1. a API calcula \(b_{2,\pm1}\) excluindo o campo proprio;
2. a forca usa o prefator dimensional correto;
3. a formula vetorial coincide com o oraculo cartesiano independente;
4. o par coincide com a solucao escalar exata `Lmax=1`;
5. os quatro valores de regressao no contato sao reproduzidos;
6. a Eq. (28) do artigo de 2026 e recuperada;
7. Silva-Bruus e recuperado nos limites apropriados;
8. todas as simetrias e escalas exigidas passam;
9. os 58 testes anteriores continuam passando;
10. toda a suite passa com warnings promovidos a erro;
11. o CSV e regenerado deterministicamente;
12. nenhuma analise de \(N\geq3\) ou extensao multipolar foi antecipada;
13. documentacao e API publica refletem corretamente o estado do projeto.

---

## 18. Verificacao final

Execute:

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python scripts/validate_t03_solver.py
python scripts/validate_t04_force.py
python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
```

Regere o CSV da T03 e confirme que ele permanece byte a byte identico ao
arquivo aprovado.

Regere duas vezes o CSV da T04 e confirme que o hash SHA-256 nao muda.

---

## 19. Relatorio, commit e push

Ao finalizar, apresente:

1. o commit-base encontrado;
2. o resultado completo de `pytest -q -W error`;
3. a contagem de testes antigos e novos;
4. o erro maximo contra o oraculo cartesiano;
5. o erro maximo contra a solucao escalar do par;
6. os quatro valores de regressao no contato;
7. o erro nos limites de recuperacao de Silva-Bruus;
8. o hash do CSV da T04;
9. a confirmacao de que o CSV da T03 nao mudou;
10. a lista de arquivos alterados;
11. `git diff --stat`;
12. a confirmacao explicita de que:
    - o solver continua em `Lmax=1`;
    - `ell=2` e somente ordem de avaliacao local;
    - nenhum resultado para \(N\geq3\) foi produzido;
    - T01 e T02 nao foram alteradas cientificamente.

Faca um commit com mensagem semelhante a:

```text
feat: compute T04 Rayleigh interaction forces
```

Depois, faca o push para a mesma branch remota atualmente utilizada.

Se o push exigir trocar o remote, sobrescrever historico, usar `--force` ou
resolver uma divergencia nao relacionada a T04, pare e informe o problema.

