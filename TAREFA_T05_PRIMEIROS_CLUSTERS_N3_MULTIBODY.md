# T05 - Primeiros clusters com \(N=3\) e contribuicao genuinamente multibody

## Instrucao de uso

Copie esta especificacao integralmente para o Codex CLI aberto na raiz do
repositorio `silva-bruus-multiscattering`.

A T05 e a primeira tarefa que deve produzir e interpretar forcas para
clusters com \(N\geq3\). O nucleo cientifico das T01-T04.1 ja foi auditado e
aprovado. Esta tarefa deve reutiliza-lo, sem alterar suas equacoes, para
comparar de forma controlada os Modelos A, B e C em trimetros.

---

## 1. Objetivo cientifico

Quantificar, pela primeira vez, a parcela da forca de interacao que nao pode
ser representada como soma de problemas isolados de duas particulas.

Para cada particula \(i\), devem ser calculadas:

\[
\mathbf F_i^A
=
\mathbf F_i^{\mathrm{SB}},
\]

\[
\mathbf F_i^B
=
\mathbf F_i^{\mathrm{2B}},
\]

\[
\mathbf F_i^C
=
\mathbf F_i^{\mathrm{MS},N=3,\ell\leq1}.
\]

A contribuicao irredutivel de tres corpos e definida por

\[
\boxed{
\Delta\mathbf F_i^{(3)}
=
\mathbf F_i^C-\mathbf F_i^B
}.
\]

Tambem deve ser registrada a correcao completa de dois corpos:

\[
\boxed{
\Delta\mathbf F_i^{(2)}
=
\mathbf F_i^B-\mathbf F_i^A
}.
\]

A decomposicao deve satisfazer, partícula por particula,

\[
\boxed{
\mathbf F_i^C-\mathbf F_i^A
=
\Delta\mathbf F_i^{(2)}
+
\Delta\mathbf F_i^{(3)}
}.
\]

Esta identidade e vetorial. Ela deve ser aplicada as componentes assinadas da
forca, e nao a modulos ou erros relativos.

O resultado central da T05 deve responder:

> Quanto da discrepancia entre Silva-Bruus e o multiple scattering Rayleigh
> global ja e explicado por correcoes completas de pares, e quanto permanece
> como efeito genuinamente coletivo de tres particulas?

---

## 2. Estado inicial aprovado

A base cientifica final da T04.1 esta na `main`, no commit:

```text
4615e7fbef3c1c4af77341d4c21edff1c967feae
```

Antes de alterar qualquer arquivo, execute:

```bash
git status --short
git branch --show-current
git rev-parse HEAD
git log -8 --oneline
```

Confirme:

1. arvore de trabalho limpa;
2. branch `main`;
3. `HEAD` contendo o commit acima;
4. nenhuma divergencia que exija sobrescrever trabalho do usuario.

Leia integralmente:

- `AGENTS.md`;
- `TASKS.md`;
- `README.md`;
- `docs/CONVENTIONS.md`;
- `docs/DECISIONS.md`;
- `docs/HANDOFF.md`;
- `TAREFA_T03_NUCLEO_MULTIPOLAR_SOLVER_RAYLEIGH.md`;
- `TAREFA_T04_FORCA_INTERACAO_MS_RAYLEIGH.md`;
- `TAREFA_T04_1_FECHAMENTO_TESTES_DOCUMENTACAO.md`;
- `src/acoustic_ms/silva_bruus.py`;
- `src/acoustic_ms/corrected_pair.py`;
- `src/acoustic_ms/solver.py`;
- `src/acoustic_ms/force.py`;
- `src/acoustic_ms/__init__.py`;
- os testes das T01-T04.1.

Prepare o ambiente e execute:

```bash
python -m pytest -q -W error
```

A linha de base obrigatoria e:

```text
85 passed
```

Confirme ainda que os CSVs aprovados continuam com os hashes:

```text
results/data/t03_solver_validation.csv
7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a

results/data/t04_pair_force_validation.csv
15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5
```

Se qualquer pre-condicao falhar, pare e relate o problema. Nao use
`git reset --hard`, `git checkout --`, `git clean`, `--force` ou outro comando
destrutivo.

---

## 3. Regime fisico e convencoes preservadas

Use exatamente o regime ja aprovado:

- fluido ideal, inviscido, homogeneo e ilimitado;
- acustica linear, harmonica e dependencia temporal \(e^{-i\omega t}\);
- esferas identicas, compressiveis, fixas e nao sobrepostas;
- onda estacionaria normalizada por

  \[
  \widetilde p^{\mathrm{ext}}=\sin(kz);
  \]

- todos os centros no plano nodal \(z=0\);
- regime de Rayleigh:

  \[
  0<ka\leq0.1;
  \]

- energia acustica explicita:

  \[
  E_0=\frac{\rho_0|v_0|^2}{4};
  \]

- vetor da fonte \(j\) para o alvo \(i\):

  \[
  \mathbf d_{ij}=\mathbf r_i-\mathbf r_j,
  \qquad
  \widehat{\mathbf d}_{ij}=\frac{\mathbf d_{ij}}{d_{ij}};
  \]

- forca radial negativa significa atracao;
- nao sobreposicao:

  \[
  d_{ij}\geq2a;
  \]

- base de harmonicos esfericos complexos ortonormais com fase de
  Condon-Shortley;
- operador de translacao orientado como `target <- source`;
- solver de espalhamento fixo em:

  \[
  L_{\max}^{\mathrm{scatter}}=1;
  \]

- avaliacao local da forca ate:

  \[
  L_{\max}^{\mathrm{evaluation}}=2.
  \]

Os modos locais \(\ell=2\) continuam sendo apenas observaveis do campo regular
reexpandido. Eles nao sao quadrupolos espalhados e nao devem realimentar o
solver.

Nao trate \(ka\), \(kd\) e \(d/a\) como parametros independentes:

\[
kd=(ka)\frac da.
\]

A T05 continua usando a forca de interacao externo-espalhado das Eqs. (22) e
(27) do artigo de 2026. Nao acrescente termos
espalhado-espalhado, forca primaria fora do plano nodal ou uma nova
integracao do tensor de tensoes.

---

## 4. Definicao operacional dos Modelos A, B e C

### 4.1 Modelo A - Silva-Bruus pairwise

Para um cluster com centros \(\{\mathbf r_i\}\),

\[
\boxed{
\mathbf F_i^A
=
\sum_{j\ne i}
\mathbf F_{i\leftarrow j}^{\mathrm{SB}}
}.
\]

Cada termo deve ser obtido pela API publica ja aprovada em
`silva_bruus.py`, preservando o sinal e a orientacao fonte-alvo.

O Modelo A contem apenas a aproximacao pairwise original de Silva-Bruus.

### 4.2 Modelo B - soma de pares isolados resolvidos por multiple scattering

Para cada par nao ordenado \(\{i,j\}\):

1. extraia as duas posicoes do cluster;
2. resolva esse par isoladamente com
   `solve_rayleigh_nodal_interaction_forces`;
3. acumule a forca retornada para as duas particulas nos indices globais
   correspondentes;
4. processe cada par exatamente uma vez.

Assim,

\[
\boxed{
\mathbf F_i^B
=
\sum_{j\ne i}
\mathbf F_{i\leftarrow j}^{C,N=2}
}.
\]

Para \(N=3\), devem ser resolvidos exatamente os tres pares:

\[
\{1,2\},\qquad\{1,3\},\qquad\{2,3\}.
\]

Cada par usa:

- o mesmo \(ka\);
- os mesmos \(f_0\) e \(f_1\);
- o mesmo \(E_0\);
- o mesmo solver \(L_{\max}=1\);
- a mesma formula de forca da T04.

### 4.3 Modelo C - multiple scattering global

Resolva as tres particulas simultaneamente por uma unica chamada:

```python
solve_rayleigh_nodal_interaction_forces(...)
```

O sistema global e

\[
(\mathbf I-\mathbf D_g\mathbf U)\mathbf s
=
\mathbf D_g\mathbf a^{\mathrm{ext}},
\]

com \(4N=12\) coeficientes complexos para o trimero e todas as cadeias de
reespalhamento permitidas por \(L_{\max}=1\).

O Modelo C admite caminhos como

\[
1\rightarrow2\rightarrow3,
\qquad
1\rightarrow3\rightarrow2\rightarrow1,
\]

que nao existem em uma soma de pares resolvidos separadamente.

---

## 5. Decisao obrigatoria: a formula da T02 nao define o Modelo B

`corrected_pair.py` implementa a expressao analitica de quinta ordem do artigo
de 2026. Ela permanece um benchmark externo importante, mas nao deve ser usada
para construir o Modelo B desta tarefa.

Se B fosse calculado pela formula truncada da T02 e C pelo solver acoplado da
T04, a diferenca

\[
\mathbf F^C-\mathbf F^B
\]

misturaria:

1. efeitos genuinamente multibody;
2. diferencas entre a expansao analitica truncada da T02 e a solucao
   numerica completa do par em \(L_{\max}=1\).

Para que \(\Delta\mathbf F^{(3)}\) tenha interpretacao unica, B e C devem usar
o mesmo solver, os mesmos coeficientes Rayleigh e o mesmo observavel de forca.
A unica diferenca deve ser:

- B: cada par e resolvido isoladamente;
- C: todas as particulas sao resolvidas simultaneamente.

Nao remova nem altere a T02. Apenas nao a use como definicao de B.

---

## 6. API cientifica a implementar

Crie um modulo de comparacao, preferencialmente:

```text
src/acoustic_ms/comparison.py
```

Implemente uma `dataclass(frozen=True)` com nome claro, por exemplo:

```python
NodalForceModelComparison
```

Ela deve expor, no minimo:

```text
model_a_forces_xy
model_b_forces_xy
model_c_forces_xy
two_body_correction_xy
irreducible_multibody_xy
global_result
```

As cinco matrizes de forca/correcao devem ter formato:

```text
(N, 2)
```

e unidade de newtons.

Implemente uma unica funcao publica principal, por exemplo:

```python
compare_nodal_force_models(
    positions_xyz,
    k,
    radius,
    energy_density,
    f0,
    f1,
    lmax=1,
)
```

Requisitos:

1. chamar o Modelo C global uma unica vez;
2. reutilizar a validacao publica ja exercida pela T04;
3. construir A apenas com a API de Silva-Bruus;
4. construir B resolvendo cada par nao ordenado uma unica vez;
5. calcular:

   ```python
   two_body_correction_xy = model_b_forces_xy - model_a_forces_xy
   irreducible_multibody_xy = model_c_forces_xy - model_b_forces_xy
   ```

6. manter o resultado completo da chamada global em `global_result`, para
   auditar residuos, condicionamento e coeficientes;
7. aceitar genericamente \(N\geq1\), embora os resultados cientificos desta
   tarefa sejam somente para \(N=3\);
8. nao introduzir cache, paralelismo ou otimizacao prematura.

Atualize `src/acoustic_ms/__init__.py` para exportar apenas as APIs publicas
necessarias.

---

## 7. Geometrias canonicas de trimetros

Crie, preferencialmente:

```text
src/acoustic_ms/geometries.py
```

As funcoes devem retornar `numpy.ndarray` com formato `(3, 3)`, coordenada
\(z=0\) e centroide na origem.

### 7.1 Cadeia linear

Para distancia de primeiros vizinhos \(d\):

\[
\mathbf r_1=(-d,0,0),\qquad
\mathbf r_2=(0,0,0),\qquad
\mathbf r_3=(d,0,0).
\]

O parametro geometrico e:

\[
d_{\min}=d.
\]

### 7.2 Triangulo equilatero

Use lado \(d\) e centroide na origem:

\[
\mathbf r_1=
\left(-\frac d2,-\frac{\sqrt3\,d}{6},0\right),
\]

\[
\mathbf r_2=
\left(\frac d2,-\frac{\sqrt3\,d}{6},0\right),
\]

\[
\mathbf r_3=
\left(0,\frac{\sqrt3\,d}{3},0\right).
\]

Aqui:

\[
d_{\min}=d.
\]

### 7.3 Triangulo escaleno de forma fixa

Antes de centralizar, use:

\[
\mathbf q_1=(0,0,0),
\qquad
\mathbf q_2=(d,0,0),
\qquad
\mathbf q_3=
\left(\frac{3d}{11},\frac{12d}{11},0\right).
\]

Subtraia o centroide:

\[
\mathbf r_i=\mathbf q_i-\frac13\sum_{n=1}^3\mathbf q_n.
\]

Os lados sao:

\[
d,\qquad
\frac{\sqrt{153}}{11}d,\qquad
\frac{\sqrt{208}}{11}d,
\]

de modo que \(d_{\min}=d\).

As funcoes de geometria devem rejeitar distancia nao finita ou nao positiva.
A restricao \(d_{\min}\geq2a\) continua pertencendo ao solver, que conhece o
raio.

Para evitar que arredondamento na construcao do triangulo equilatero transforme
contato exato em uma sobreposicao espuria, os resultados de producao da T05
devem iniciar em:

\[
\frac{d_{\min}}a=2.1,
\]

e nao exatamente em \(2\).

Nao altere a tolerancia de nao sobreposicao do solver nesta tarefa.

---

## 8. Oraculo independente reduzido para o plano nodal

Os testes de \(N=3\) nao podem validar o codigo chamando apenas as mesmas
matrizes de traducao e a mesma formula multipolar da producao.

No plano nodal, dentro de \(L_{\max}=1\), apenas \(s_{10}^{(i)}\) sobrevive.
Construa em `tests/test_multibody.py` um oraculo exclusivo de teste que resolva
diretamente um sistema escalar \(N\times N\).

Defina:

\[
s_1=i\frac{f_1}{6}(ka)^3,
\qquad
a_{10}^{\mathrm{ext}}=\sqrt{12\pi}.
\]

Para \(i\ne j\):

\[
G_{ij}
=
h_0^{(1)}(kd_{ij})
+
h_2^{(1)}(kd_{ij}).
\]

Monte:

\[
M_{ij}
=
\delta_{ij}
-
(1-\delta_{ij})s_1G_{ij},
\]

e resolva:

\[
\boxed{
\mathbf M\boldsymbol{\sigma}
=
s_1a_{10}^{\mathrm{ext}}\mathbf 1
}.
\]

Cada componente \(\sigma_i\) e a referencia independente para
\(s_{10}^{(i)}\).

Para a forca, defina:

\[
Q(x)
=
\frac{h_1^{(1)\prime}(x)}{x}
-
\frac{h_1^{(1)}(x)}{x^2}.
\]

Entao:

\[
\boxed{
\mathbf F_i^{C,\mathrm{oracle}}
=
4\pi k a^3E_0
\sqrt{\frac{3}{4\pi}}
\sum_{j\ne i}
\operatorname{Re}
\left[
f_1^*\sigma_jQ(kd_{ij})
\right]
\widehat{\mathbf d}_{ij,xy}
}.
\]

Regras de independencia:

- use diretamente `scipy.special.spherical_jn` e
  `scipy.special.spherical_yn`;
- nao chame `translation_matrix`;
- nao chame `spherical_hankel1` do pacote;
- nao chame `solve_rayleigh_nodal`;
- nao chame `solve_rayleigh_nodal_interaction_forces`;
- nao reutilize \(b_{2,\pm1}\);
- o oraculo deve permanecer apenas nos testes.

Compare:

1. \(\sigma_i\) com os coeficientes \(s_{10}^{(i)}\) do solver global;
2. a forca cartesiana do oraculo com `model_c_forces_xy`;
3. os demais modos de espalhamento do solver com zero numerico.

Use tolerancias da ordem de:

```text
rtol = 3e-12
atol = 3e-13
```

ajustando apenas se houver justificativa numerica registrada.

---

## 9. Metricas

Crie, preferencialmente:

```text
src/acoustic_ms/metrics.py
```

### 9.1 Erro simetrico por particula

Para uma referencia \(\mathbf F_i^{\mathrm{ref}}\) e um modelo
\(\mathbf F_i^{\mathrm{mod}}\):

\[
\varepsilon_{i,\mathrm{sym}}
=
\frac{
2\left\|
\mathbf F_i^{\mathrm{ref}}-\mathbf F_i^{\mathrm{mod}}
\right\|
}{
\left\|\mathbf F_i^{\mathrm{ref}}\right\|
+
\left\|\mathbf F_i^{\mathrm{mod}}\right\|
}.
\]

Use C como referencia e calcule a metrica para A e B.

A metrica deve permanecer em:

\[
0\leq\varepsilon_{i,\mathrm{sym}}\leq2.
\]

Quando ambos os vetores forem numericamente nulos, retorne zero. A deteccao
de zero deve ser relativa a escala global das forcas da configuracao, usando
um limiar da ordem de algumas dezenas ou centenas de
`numpy.finfo(float).eps`; nao trate residuos de \(10^{-16}\) em uma
configuracao de escala \(O(1)\) como uma mudanca de direcao fisica.

### 9.2 Erro RMS global

\[
\varepsilon_{\mathrm{RMS}}
=
\left[
\frac{
\sum_i
\left\|
\mathbf F_i^{\mathrm{ref}}-\mathbf F_i^{\mathrm{mod}}
\right\|^2
}{
\sum_i
\left\|\mathbf F_i^{\mathrm{ref}}\right\|^2
}
\right]^{1/2}.
\]

Novamente, use C como referencia e avalie A e B.

Se referencia e modelo forem ambos nulos, retorne zero. Se a referencia for
nula mas o modelo nao for, retorne infinito ou lance uma excecao
documentada; nao permita `0/0` silencioso.

### 9.3 Erro angular

Implemente uma metrica angular por particula:

\[
\theta_i
=
\arccos
\left[
\frac{
\mathbf F_i^{\mathrm{ref}}\cdot\mathbf F_i^{\mathrm{mod}}
}{
\|\mathbf F_i^{\mathrm{ref}}\|
\|\mathbf F_i^{\mathrm{mod}}\|
}
\right].
\]

Retorne o resultado explicitamente em graus e use `clip` no argumento do
`arccos`. Quando um dos vetores for numericamente nulo, retorne `NaN`, pois a
direcao nao esta definida.

Teste as metricas com vetores identicos, opostos, ortogonais, ambos nulos e
com um unico vetor nulo.

---

## 10. Testes cientificos obrigatorios

Crie:

```text
tests/test_multibody.py
```

Nao duplique testes antigos sem necessidade.

### 10.1 Reducoes \(N=1\) e \(N=2\)

Para \(N=1\):

\[
\mathbf F^A=\mathbf F^B=\mathbf F^C=\mathbf0.
\]

Para \(N=2\):

\[
\boxed{
\mathbf F^B=\mathbf F^C
}
\]

dentro da precisao numerica, pois B possui exatamente um par e deve chamar a
mesma API da T04.

Consequentemente:

\[
\Delta\mathbf F^{(3)}=\mathbf0
\qquad (N\leq2).
\]

### 10.2 Oraculo escalar independente

Para as tres geometrias canonicas:

- compare todos os \(s_{10}^{(i)}\);
- confirme que os outros tres modos do solver \(L_{\max}=1\) sao nulos;
- compare todas as componentes de \(\mathbf F_i^C\);
- cubra ao menos um caso obliquo/escaleno sem simetria axial.

### 10.3 Identidade telescopica

Verifique elemento a elemento:

\[
\mathbf F^C-\mathbf F^A
=
(\mathbf F^B-\mathbf F^A)
+
(\mathbf F^C-\mathbf F^B).
\]

### 10.4 Permutacao das particulas

Use um trimero escaleno e teste mais de uma permutacao, incluindo uma
permutacao nao trivial como:

```python
order = np.array([2, 0, 1])
```

Depois de aplicar a mesma permutacao aos resultados de referencia, compare:

- A;
- B;
- C;
- \(\Delta\mathbf F^{(2)}\);
- \(\Delta\mathbf F^{(3)}\);
- coeficientes do solver global.

### 10.5 Translacao e rotacao

Para um trimero escaleno:

- translade todos os centros pelo mesmo vetor no plano;
- rotacione todos os centros por um angulo que nao seja multiplo de
  \(\pi/2\);
- A, B, C e as duas correcoes devem ser invariantes por translacao;
- os vetores devem transformar covariantemente sob rotacao.

### 10.6 Escalamento dimensional

Mantenha \(ka\) e todas as razoes \(d_{ij}/a\) fixas:

\[
a\rightarrow\lambda a,
\qquad
\mathbf r_i\rightarrow\lambda\mathbf r_i,
\qquad
k\rightarrow\frac{k}{\lambda}.
\]

Com \(E_0\) fixo:

\[
\mathbf F\rightarrow\lambda^2\mathbf F.
\]

Verifique isso para A, B, C, \(\Delta\mathbf F^{(2)}\) e
\(\Delta\mathbf F^{(3)}\).

### 10.7 Energia, contraste e plano nodal

Verifique:

- linearidade em \(E_0\);
- forca nula para \(E_0=0\);
- forca nula para \(f_1=0\);
- independencia de \(f_0\) para trimetros no plano nodal;
- rejeicao das entradas invalidas ja cobertas pela API da T04;
- rejeicao de `lmax != 1`.

### 10.8 Simetrias geometricas

#### Cadeia linear

- a particula central deve ter forca numericamente nula;
- as particulas externas devem ter forcas opostas;
- nao deve haver componente transversal fisica;
- as correcoes de dois e tres corpos devem respeitar a mesma simetria.

#### Triangulo equilatero

- as tres forcas devem ter o mesmo modulo;
- cada forca deve estar na direcao radial relativa ao centroide;
- a configuracao deve ser covariante por rotacao de \(120^\circ\);
- a soma vetorial deve ser numericamente nula por simetria.

#### Triangulo escaleno

Nao imponha artificialmente:

\[
\sum_i\mathbf F_i^C=\mathbf0.
\]

Os Modelos A e B sao somas de pares antissimetricos e, portanto, possuem soma
nula ate arredondamento. Essa propriedade nao decorre automaticamente para o
Modelo C, pois ele e um observavel de interacao externo-espalhado em um sistema
coletivo alimentado pelo campo acustico externo. Para o escaleno, valide a
forca e sua soma pelo oraculo independente, e nao por uma hipotese de
acao-reacao par a par. Uma soma nao nula neste observavel nao deve ser chamada
automaticamente de forca total ou de recuo total do cluster, pois a T05 nao
calcula a forca irrestrita com todos os termos de momento.

### 10.9 Limite de acoplamento fraco

Para separacoes crescentes ou \(f_1\) decrescente, verifique numericamente:

\[
\mathbf F^B-\mathbf F^A\rightarrow\mathbf0,
\]

\[
\mathbf F^C-\mathbf F^B\rightarrow\mathbf0.
\]

Nao imponha monotonicidade ponto a ponto em uma faixa que possa conter
oscilacoes de \(\cos(kd)\). Use pontos selecionados no regime em que a
comparacao seja bem condicionada e documente a metrica usada.

---

## 11. Regressoes numericas canonicas

Use:

\[
a=1,\qquad
E_0=1,\qquad
ka=0.1,\qquad
f_0=0,\qquad
f_1=0.8.
\]

Os valores numericos correspondem simultaneamente a
\(\mathbf F/(a^2E_0)\), pois \(a^2E_0=1\).

### 11.1 Resumo das tres geometrias

| Geometria | \(d_{\min}/a\) | \(\varepsilon_{\mathrm{RMS}}(A,C)\) | \(\varepsilon_{\mathrm{RMS}}(B,C)\) |
|---|---:|---:|---:|
| cadeia linear | 2.1 | 0.0831973786371517 | 0.0430867438601629 |
| equilatero | 2.1 | 0.0882674934567006 | 0.0461714463617981 |
| escaleno fixo | 2.2 | 0.0578621870496697 | 0.0270655878714466 |

### 11.2 Cadeia linear

Para as posicoes:

\[
(-2.1,0,0),\quad(0,0,0),\quad(2.1,0,0),
\]

o Modelo C deve produzir, desprezando ruido transversal:

\[
\mathbf F^C_1
=
(0.725204866335105,\ 0),
\]

\[
\mathbf F^C_2=\mathbf0,
\]

\[
\mathbf F^C_3
=
(-0.725204866335105,\ 0).
\]

Para a particula externa esquerda:

\[
F^A_{1x}=0.664869722481118,
\]

\[
F^B_{1x}=0.693958150013181,
\]

\[
\Delta F^{(3)}_{1x}=0.0312467163219243.
\]

### 11.3 Triangulo equilatero

Para \(d/a=2.1\), a particula superior deve ter:

\[
\mathbf F^C_3
\simeq
(0,-1.187167981785305).
\]

As outras duas forcas devem ser as rotacoes por \(\pm120^\circ\):

\[
\mathbf F^C_1
\simeq
(1.028117630785575,\ 0.593583990892652),
\]

\[
\mathbf F^C_2
\simeq
(-1.028117630785576,\ 0.593583990892652).
\]

### 11.4 Triangulo escaleno

Use, antes da centralizacao:

\[
(0,0,0),\qquad(2.2,0,0),\qquad(0.6,2.4,0).
\]

A centralizacao nao altera as forcas. O Modelo C deve produzir:

\[
\mathbf F^C
\simeq
\begin{pmatrix}
0.632645020072720 & 0.330698096562575\\
-0.657831717956925 & 0.154039050571838\\
0.019458368752579 & -0.493350707384327
\end{pmatrix}.
\]

A soma resultante e:

\[
\sum_i\mathbf F_i^C
\simeq
(-0.005728329131626,\ -0.008613560249914).
\]

Esse valor nao deve ser forcado a zero. Ele deve concordar com o oraculo
escalar independente.

Para essa configuracao:

\[
\Delta\mathbf F^{(3)}
\simeq
\begin{pmatrix}
0.011556768765901 & 0.006155893167571\\
-0.018110136050905 & 0.004385728399909\\
0.000825038153378 & -0.019155181817394
\end{pmatrix}.
\]

Use tolerancias relativas da ordem de \(5\times10^{-12}\) e absolutas
compativeis com a escala. Nao transforme ruido numerico \(O(10^{-16})\) em
regressao fisica.

---

## 12. Varredura de producao da T05

Crie:

```text
scripts/validate_t05_trimers.py
```

O script deve usar apenas APIs publicas do pacote para A, B e C.

### 12.1 Parametros

Fixe:

\[
a=1,\qquad
E_0=1,\qquad
ka=0.1,\qquad
f_0=0.
\]

Use:

\[
f_1\in\{0.1,\ 0.4,\ 0.8,\ 1.0\}.
\]

Para cada uma das tres familias geometricas, use exatamente:

\[
\frac{d_{\min}}a
\in[2.1,10.0],
\]

com:

```python
np.linspace(2.1, 10.0, 160)
```

Use o mesmo vetor de 160 separacoes nas tres geometrias.

### 12.2 Quantidades por configuracao

Registre, no minimo:

- geometria;
- \(ka\);
- \(f_0\);
- \(f_1\);
- \(d_{\min}/a\);
- \(\varepsilon_{\mathrm{RMS}}(A,C)\);
- \(\varepsilon_{\mathrm{RMS}}(B,C)\);
- maximo do erro simetrico por particula para A contra C;
- maximo do erro simetrico por particula para B contra C;
- maximo erro angular definido para A contra C;
- maximo erro angular definido para B contra C;
- norma RMS de \(\Delta\mathbf F^{(2)}/(a^2E_0)\);
- norma RMS de \(\Delta\mathbf F^{(3)}/(a^2E_0)\);
- resíduo relativo do solver global;
- numero de condicao do sistema global;
- componentes da soma \(\sum_i\mathbf F_i^C/(a^2E_0)\).

Para angulos indefinidos, grave `NaN` de forma deterministica e documentada.

### 12.3 Arquivos

Gere:

```text
results/data/t05_trimer_regression.csv
results/data/t05_trimer_sweep.csv
results/figures/t05_trimer_model_errors.png
```

`t05_trimer_regression.csv` deve conter, por partícula, os valores A, B, C,
\(\Delta\mathbf F^{(2)}\) e \(\Delta\mathbf F^{(3)}\) para os tres casos
canonicos da Secao 11.

`t05_trimer_sweep.csv` deve conter a varredura completa da Secao 12.

### 12.4 Figura

Produza uma figura com tres linhas e duas colunas:

- uma linha para cadeia, equilatero e escaleno;
- coluna esquerda:

  \[
  100\,\varepsilon_{\mathrm{RMS}}(A,C);
  \]

- coluna direita:

  \[
  100\,\varepsilon_{\mathrm{RMS}}(B,C).
  \]

Em cada painel, mostre as quatro curvas de \(f_1\). Use:

- eixo horizontal \(d_{\min}/a\);
- eixo vertical em porcentagem;
- escala logaritmica em \(y\) se todos os valores plotados forem positivos;
- mesmas cores para o mesmo \(f_1\) em todos os paineis;
- rotulos, unidades adimensionais, legenda e grade discreta;
- titulo ou legenda que identifique A como Silva-Bruus, B como soma de pares
  MS e C como MS global.

A figura deve comunicar claramente que:

- A contra C contem correcoes de dois corpos e multibody;
- B contra C isola o efeito multibody no truncamento Rayleigh.

Nao afirme convergencia multipolar ou validade do Modelo D.

---

## 13. Determinismo e auditoria dos artefatos

Execute o script duas vezes.

Depois da primeira execucao, registre os hashes SHA-256 dos dois CSVs e do
PNG. Depois da segunda, confirme que os hashes sao identicos.

Os CSVs devem:

- usar cabecalhos explicitos;
- ter ordem deterministica de geometrias, contrastes, separacoes e
  particulas;
- usar precisao suficiente para regressao;
- nao depender de estado aleatorio;
- nao conter indices de DataFrame;
- terminar com newline.

Verifique programaticamente:

- nenhuma forca ou metrica obrigatoria infinita nas configuracoes de
  producao;
- `NaN` apenas nos erros angulares matematicamente indefinidos;
- residuos finitos;
- numeros de condicao finitos;
- identidade telescopica em todas as linhas de regressao;
- \(0\leq\varepsilon_{\mathrm{sym}}\leq2\).

---

## 14. Arquivos permitidos

E permitido criar:

```text
src/acoustic_ms/comparison.py
src/acoustic_ms/geometries.py
src/acoustic_ms/metrics.py
tests/test_multibody.py
scripts/validate_t05_trimers.py
results/data/t05_trimer_regression.csv
results/data/t05_trimer_sweep.csv
results/figures/t05_trimer_model_errors.png
TAREFA_T05_PRIMEIROS_CLUSTERS_N3_MULTIBODY.md
```

E permitido atualizar:

```text
src/acoustic_ms/__init__.py
README.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
```

Nao altere:

```text
src/acoustic_ms/force.py
src/acoustic_ms/solver.py
src/acoustic_ms/translation.py
src/acoustic_ms/scattering.py
src/acoustic_ms/special.py
src/acoustic_ms/gaunt.py
src/acoustic_ms/silva_bruus.py
src/acoustic_ms/corrected_pair.py
tests/test_force.py
```

Nao altere os CSVs e a figura das tarefas anteriores.

Se um novo teste revelar um defeito real em um modulo cientifico protegido,
nao o corrija silenciosamente. Pare e apresente:

1. o menor contraexemplo;
2. a saida numerica;
3. o modulo afetado;
4. o impacto sobre T01-T04.1;
5. a menor correcao proposta.

---

## 15. Documentacao obrigatoria

### `docs/CONVENTIONS.md`

Documente:

- as definicoes exatas de A, B e C;
- \(\Delta\mathbf F^{(2)}\) e
  \(\Delta\mathbf F^{(3)}\);
- o uso de \(\mathbf F/(a^2E_0)\) para resultados normalizados;
- as definicoes das metricas e o tratamento de forcas numericamente nulas;
- que B usa pares resolvidos pela T04, e nao a formula truncada da T02.

### `docs/DECISIONS.md`

Registre as decisoes metodologicas:

- T05 produz apenas resultados \(N=3\);
- B e C compartilham o mesmo solver/observavel;
- a diferenca C-B isola multibody dentro de \(L_{\max}=1\);
- o oraculo escalar e exclusivo de teste;
- nao se impoe soma de forcas nula para o escaleno global;
- a T05 nao mede correcao multipolar.

### `docs/HANDOFF.md`

Registre:

- arquivos alterados;
- equacoes implementadas;
- parametros das tres geometrias;
- contagem real de testes;
- comandos executados;
- erros maximos contra o oraculo;
- regressões da Secao 11;
- residuos e condicionamento maximos da varredura;
- hashes dos tres novos artefatos;
- confirmacao de hashes preservados das T03 e T04;
- limitacoes e proximo passo.

### `TASKS.md`

Adicione e marque como concluida somente depois de todos os criterios:

```text
- [x] T05 -- first N=3 Model A/B/C comparison and irreducible multibody force.
```

### `README.md`

Atualize o resumo do estado do projeto para nao afirmar que multiple
scattering ou forcas \(N\geq3\) ainda estao fora do escopo. Descreva com
precisao que o estado atual, apos T05, chega a trimetros em
\(L_{\max}=1\), sem Modelo D.

---

## 16. Verificacoes finais

Execute, no minimo:

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python -m pytest -q -W error
python scripts/validate_t03_solver.py
python scripts/validate_t04_force.py
python scripts/validate_t05_trimers.py
python scripts/validate_t05_trimers.py
git diff --check
git status --short
git diff --stat
```

Calcule:

```bash
sha256sum \
  results/data/t03_solver_validation.csv \
  results/data/t04_pair_force_validation.csv \
  results/data/t05_trimer_regression.csv \
  results/data/t05_trimer_sweep.csv \
  results/figures/t05_trimer_model_errors.png
```

Inspecione visualmente o PNG final. Verifique:

- textos legiveis;
- eixos nao cortados;
- legenda sem sobreposicao;
- curvas distinguiveis;
- ausencia de paineis vazios;
- ausencia de `NaN` ou `inf` renderizado como dado fisico.

Nao declare antecipadamente a contagem final de testes. Registre o total
efetivamente coletado.

---

## 17. Criterios de aprovacao

A T05 somente estara concluida se:

1. os 85 testes anteriores continuarem passando;
2. os novos testes passarem com warnings tratados como erros;
3. A, B e C estiverem definidos e implementados sem ambiguidade;
4. B usar pares isolados resolvidos pela T04;
5. C usar uma unica solucao global por configuracao;
6. a identidade telescopica for satisfeita;
7. B e C coincidirem para \(N=2\);
8. o oraculo escalar independente concordar com coeficientes e forcas de C
   para as tres geometrias;
9. permutacao, translacao, rotacao e escalamento forem verificados;
10. as simetrias da cadeia e do equilatero forem satisfeitas;
11. o escaleno concordar com o oraculo sem impor soma de forcas nula;
12. as regressoes numericas da Secao 11 forem reproduzidas;
13. os dois CSVs e o PNG forem deterministicos;
14. os hashes das T03 e T04 permanecerem inalterados;
15. nenhum modulo cientifico protegido for modificado;
16. a documentacao registrar claramente o que C-B mede;
17. nenhum resultado de \(N>3\) ou \(L_{\max}>1\) for produzido.

---

## 18. Fora do escopo

Nao implemente nesta tarefa:

- Modelo D;
- \(L_{\max}=2,3\) ou \(5\) no espalhamento;
- coeficientes exatos de esfera fora da aproximacao Rayleigh;
- clusters com \(N=4,6,10,\ldots\);
- aneis, redes, clusters aleatorios ou geometrias tridimensionais;
- trajetorias ou dinamica de particulas;
- torque;
- viscosidade;
- streaming acustico;
- paredes;
- interacao hidrodinamica;
- termos espalhado-espalhado;
- raio espectral \(\rho(\mathbf K)\);
- criterio \(\Lambda_i\);
- GMRES, `LinearOperator`, Numba, Joblib, GPU ou FMM;
- notebooks como implementacao cientifica principal;
- ajuste empirico de um criterio de validade.

Esses itens pertencem a tarefas posteriores.

---

## 19. Commit, push e relatorio de retorno

Antes do commit:

```bash
git diff --check
git status --short
git diff --stat
```

Confirme que nenhum arquivo fora do escopo foi alterado.

Use a mensagem:

```text
feat: compare three-particle multibody forces
```

Faca o push para a mesma `main` remota ja utilizada.

Ao finalizar, apresente:

```markdown
# Relatorio T05

## Commit e push

## Arquivos criados ou alterados

## Definicao implementada dos Modelos A, B e C

## Oraculo independente

## Testes executados e contagem final

## Regressoes numericas

## Resultados das tres geometrias

## Residuos e condicionamento

## Hashes dos artefatos

## Confirmacao de preservacao das T01-T04.1

## Limitacoes

## Proxima acao sugerida
```

Inclua o hash completo do commit, a saida de
`python -m pytest -q -W error`, os hashes dos artefatos e um resumo objetivo
do `git diff --stat`.

Nao avance automaticamente para \(N>3\), Modelo D ou criterios de validade.
