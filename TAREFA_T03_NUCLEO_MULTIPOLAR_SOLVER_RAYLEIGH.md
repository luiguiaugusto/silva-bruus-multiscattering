# T03 - Núcleo multipolar e solver acoplado Rayleigh (\(L_{\max}=1\))

## Instrução de uso

Copie esta especificação integralmente para o Codex CLI aberto na raiz do
repositório `silva-bruus-multiscattering`.

Esta tarefa parte da `main` aprovada após a T02. Antes de alterar qualquer
arquivo:

1. execute `git status --short` e preserve alterações preexistentes;
2. leia `AGENTS.md`, `README.md`, `TASKS.md`, `docs/CONVENTIONS.md`,
   `docs/DECISIONS.md` e `docs/HANDOFF.md`;
3. leia todo o código e todos os testes existentes;
4. consulte no artigo de 2026 as Eqs. (7)-(16), (23)-(24), (A1) e os
   Apêndices B-C;
5. não altere as fórmulas nem os resultados aprovados nas T01 e T02.

Não faça `commit` nem `push`, a menos que o usuário solicite isso
separadamente.

---

## 1. Objetivo científico

Implemente o primeiro solver de múltiplo espalhamento do projeto para
partículas no plano nodal, no regime de Rayleigh e com

\[
\boxed{L_{\max}=1.}
\]

Cada partícula será representada pelos modos

\[
(\ell,m)=(0,0),(1,-1),(1,0),(1,1),
\]

totalizando quatro coeficientes por partícula. Para \(N_p\) partículas, o
sistema completo terá \(4N_p\) incógnitas complexas.

O solver deverá obter os coeficientes espalhados acoplados

\[
s_{nm}^{(i)}
=s_n
\left[
a_{nm}^{(i)}
+\sum_{j\ne i}\sum_{n'm'}
\mathcal S_{nm,n'm'}^{\,j\to i}s_{n'm'}^{(j)}
\right]
\]

por meio de uma solução linear. A inversão do sistema deve ressomar
implicitamente todas as cadeias de reespalhamento permitidas pela base
\(\ell\leq1\).

### Resultado da T03

A T03 entrega:

- base e indexação multipolar;
- funções esféricas e harmônicos esféricos;
- coeficientes de Gaunt;
- operador de translação entre centros;
- coeficientes de espalhamento Rayleigh de uma esfera;
- coeficientes da onda estacionária no plano nodal;
- montagem e solução do sistema acoplado para \(N_p\) partículas;
- benchmarks independentes para uma e duas partículas.

### O que a T03 não entrega

A T03 **não** deve:

- calcular força de radiação;
- comparar forças com Silva-Bruus ou com a T02;
- produzir ainda a contribuição irredutível de três corpos;
- implementar trajetórias ou dinâmica de partículas;
- avançar o solver de produção para \(L_{\max}=3\) ou \(5\);
- incluir paredes, viscosidade, streaming ou interação hidrodinâmica.

Sem uma fórmula de força acoplada, os coeficientes \(s_{nm}^{(i)}\) não devem
ser apresentados como uma nova previsão de força.

---

## 2. Escopo físico

Considere:

- fluido ideal, homogêneo, invíscido e ilimitado;
- acústica linear, harmônica e tridimensional;
- convenção temporal

  \[
  p(\mathbf r,t)=\operatorname{Re}
  \{p(\mathbf r)e^{-i\omega t}\};
  \]

- funções de onda outgoing representadas por \(h_\ell^{(1)}\);
- esferas compressíveis, idênticas e de raio \(a\);
- centros fixos no plano nodal \(z=0\);
- onda estacionária externa

  \[
  \widetilde p^{\mathrm{ext}}(\mathbf r)=\sin(kz);
  \]

- regime de Rayleigh

  \[
  0<ka\leq0.1;
  \]

- não sobreposição

  \[
  r_{ij}\geq2a.
  \]

O solver trabalha com o campo de pressão normalizado usado no artigo de 2026.
Nesta tarefa não entra \(E_0\), pois nenhuma força será calculada.

Os contrastes continuam sendo

\[
f_0=1-\frac{\kappa_p}{\kappa_0},
\qquad
f_1=\frac{2(\rho_p/\rho_0-1)}
{2\rho_p/\rho_0+1}.
\]

Aceite \(f_0\) e \(f_1\) como argumentos explícitos. Valide \(f_1\) no
intervalo físico

\[
-2\leq f_1\leq1.
\]

Não imponha um limite inferior artificial a \(f_0\), pois partículas muito
compressíveis podem apresentar \(f_0\ll-1\). Exija apenas que \(f_0\) seja real
e finito.

---

## 3. Convenções matemáticas obrigatórias

### 3.1 Harmônicos esféricos

Use harmônicos esféricos complexos ortonormais com fase de Condon-Shortley:

\[
Y_\ell^m(\theta,\phi)
=
\sqrt{
\frac{2\ell+1}{4\pi}
\frac{(\ell-m)!}{(\ell+m)!}
}
P_\ell^m(\cos\theta)e^{im\phi}.
\]

Os ângulos são:

- \(\theta\): ângulo polar ou colatitude, \(0\leq\theta\leq\pi\);
- \(\phi\): azimute, \(0\leq\phi<2\pi\).

Deve valer

\[
Y_\ell^{-m}=(-1)^mY_\ell^{m*}.
\]

Use `scipy.special.sph_harm_y`, cuja ordem de argumentos é
`(ell, m, theta, phi)`. Não use uma API antiga com os ângulos trocados.

### 3.2 Ondas radiais

Defina

\[
h_\ell^{(1)}(x)=j_\ell(x)+iy_\ell(x).
\]

Com \(e^{-i\omega t}\), essa é a solução outgoing. Implemente também as
derivadas em relação ao argumento \(x\), usando as rotinas da SciPy.

### 3.3 Ordem dos modos

Use ordenação por partícula e, dentro de cada partícula, por \(\ell\) crescente
e \(m=-\ell,\ldots,\ell\):

\[
\operatorname{mode\_index}(\ell,m)=\ell^2+\ell+m.
\]

Assim, para \(L_{\max}=1\):

| índice | \((\ell,m)\) |
|---:|:---:|
| 0 | \((0,0)\) |
| 1 | \((1,-1)\) |
| 2 | \((1,0)\) |
| 3 | \((1,1)\) |

O índice global é

\[
\alpha=i(L_{\max}+1)^2+\operatorname{mode\_index}(\ell,m),
\]

com partículas indexadas internamente a partir de zero.

Documente essa convenção em `docs/CONVENTIONS.md`.

---

## 4. Funções especiais e coeficientes de Gaunt

Implemente o coeficiente

\[
\mathcal G(n',m';q,\mu;n,m)
=
\int
Y_{n'}^{m'}(\Omega)
Y_q^\mu(\Omega)
Y_n^{-m}(\Omega)\,d\Omega,
\]

com

\[
\mu=m-m'.
\]

Use a relação com símbolos \(3j\):

\[
\begin{aligned}
\mathcal G
={}&
\sqrt{
\frac{(2n'+1)(2q+1)(2n+1)}{4\pi}
}
\\
&\times
\begin{pmatrix}
n'&q&n\\
0&0&0
\end{pmatrix}
\begin{pmatrix}
n'&q&n\\
m'&\mu&-m
\end{pmatrix}.
\end{aligned}
\]

É permitido usar `sympy.physics.wigner.gaunt` ou `wigner_3j` para gerar o valor
e convertê-lo para ponto flutuante. Armazene resultados com cache, pois os
mesmos acoplamentos serão reutilizados.

Retorne zero imediatamente quando alguma regra de seleção falhar:

\[
m'+\mu-m=0,
\]

\[
|n'-q|\leq n\leq n'+q,
\]

\[
n'+q+n\ \text{par},
\]

além dos limites usuais dos números azimutais.

Não use integração numérica para construir o operador de produção. A
quadratura será usada somente como teste independente.

---

## 5. Operador de translação

### 5.1 Definição operacional do deslocamento

Para traduzir o campo outgoing da partícula-fonte \(j\) para uma expansão
regular no centro da partícula-alvo \(i\), defina no código

\[
\boxed{
\mathbf R_{i\leftarrow j}=\mathbf r_j-\mathbf r_i
}
\]

isto é, o vetor que aponta do centro-alvo para o centro-fonte.

Essa definição deve ficar explícita na API. A função pública deve receber
`source_position` e `target_position`, e não um vetor cujo sentido fique
ambíguo.

### 5.2 Coeficiente de separação

As linhas correspondem ao modo-alvo \((n,m)\) e as colunas ao modo-fonte
\((n',m')\). Implemente

\[
\boxed{
\mathcal S_{nm,n'm'}^{\,j\to i}
=
4\pi(-1)^{n+n'+m}i^{\,n-n'}
\sum_{\ell_q=0}^{(n+n'-q_0)/2}
i^q h_q^{(1)}(kR_{i\leftarrow j})
\left[
Y_q^{m-m'}(\Omega_{i\leftarrow j})
\right]^*
\mathcal G(n',m';q,m-m';n,m)
}
\]

com

\[
q=q_0+2\ell_q,
\]

e \(q_0\) sendo o menor inteiro que satisfaz

\[
q_0\geq\max(|n-n'|,|m-m'|),
\qquad
q_0\bmod2=(n+n')\bmod2.
\]

O conjugado do harmônico e o fator
\((-1)^{n+n'+m}\) nesta forma são necessários porque a implementação adota
harmônicos **complexos** com fase de Condon-Shortley e definiu
\(\mathbf R_{i\leftarrow j}\) no sentido alvo-para-fonte. Eles convertem a
forma da Eq. (13) para essas convenções computacionais. Não remova esses
fatores por comparação visual isolada com a equação impressa; a autoridade
numérica final é o teste direto do teorema de reexpansão da Seção 5.3.

A matriz deve obedecer

\[
\mathbf b^{(j\to i)}
=
\mathbf U_{ij}\mathbf s^{(j)},
\]

onde \(\mathbf b^{(j\to i)}\) contém coeficientes regulares no centro \(i\).

### 5.3 Teste físico obrigatório do teorema de reexpansão

Não aceite o operador apenas com testes internos de Gaunt ou de simetria.
Compare diretamente:

\[
h_{n'}^{(1)}(k|\mathbf r-\mathbf r_j|)
Y_{n'}^{m'}(\Omega_j)
\]

com

\[
\sum_{n=0}^{L_{\mathrm{teste}}}
\sum_{m=-n}^{n}
\mathcal S_{nm,n'm'}^{\,j\to i}
j_n(k|\mathbf r-\mathbf r_i|)
Y_n^m(\Omega_i),
\]

para:

- pelo menos dois modos-fonte, incluindo um com \(m'\neq0\);
- deslocamento entre centros não alinhado com o eixo \(z\);
- pontos de avaliação com

  \[
  |\mathbf r-\mathbf r_i|<|\mathbf r_j-\mathbf r_i|;
  \]

- pelo menos um ponto com \(x\), \(y\) e \(z\) todos não nulos;
- \(L_{\mathrm{teste}}\geq10\).

A meta é erro relativo menor que \(10^{-9}\) para pontos suficientemente
próximos do centro-alvo. Esse teste é obrigatório para detectar troca do
sentido do deslocamento, dos ângulos, de conjugação ou das fases de \(i\).

Como teste analítico adicional, para centros separados no plano \(xy\),

\[
\boxed{
\mathcal S_{10,10}
=h_0^{(1)}(kd)+h_2^{(1)}(kd).
}
\]

---

## 6. Espalhamento Rayleigh de uma esfera

Nesta tarefa, use apenas os termos dominantes de monopolo e dipolo:

\[
\boxed{
s_0=-i\frac{f_0}{3}(ka)^3
}
\]

e

\[
\boxed{
s_1=i\frac{f_1}{6}(ka)^3.
}
\]

A matriz de espalhamento de cada esfera é diagonal:

\[
\mathbf D
=
\operatorname{diag}(s_0,s_1,s_1,s_1).
\]

Não implemente nesta tarefa os coeficientes \(s_3\) e \(s_5\) do Apêndice A.
Eles pertencem à futura extensão multipolar.

Teste, em particular, que

\[
f_0=0\Rightarrow s_0=0,
\qquad
f_1=0\Rightarrow s_1=0,
\]

e que os coeficientes apresentam escala \((ka)^3\).

---

## 7. Campo externo no plano nodal

Para

\[
\widetilde p^{\mathrm{ext}}=\sin(kz),
\]

implemente os beam-shape coefficients da Eq. (23) do artigo de 2026:

\[
a_{n0}=
\begin{cases}
0,&n\ \text{par},\\[2mm]
(-1)^{(n-1)/2}\sqrt{4\pi(2n+1)},&n\ \text{ímpar},
\end{cases}
\]

e

\[
a_{nm}=0,\qquad m\neq0.
\]

Para \(L_{\max}=1\):

\[
\boxed{
a_{10}=\sqrt{12\pi}
}
\]

e todos os demais coeficientes externos são nulos.

Implemente a função de maneira geral para qualquer \(L_{\max}\), mas o solver
da T03 deve rejeitar \(L_{\max}\neq1\).

---

## 8. Sistema global acoplado

Para cada par \(i\neq j\), construa o bloco

\[
\mathbf U_{ij},
\]

com linhas associadas aos modos-alvo em \(i\) e colunas associadas aos
modos-fonte em \(j\). Use

\[
\mathbf U_{ii}=\mathbf0.
\]

Monte a matriz de espalhamento global em blocos:

\[
\mathbf D_g
=
\operatorname{blockdiag}(\mathbf D,\ldots,\mathbf D).
\]

Empilhando os coeficientes por partícula, a Eq. (16) deve ser escrita como

\[
\mathbf s
=
\mathbf D_g
\left(
\mathbf a^{\mathrm{ext}}+\mathbf U\mathbf s
\right).
\]

Portanto,

\[
\boxed{
\left(\mathbf I-\mathbf D_g\mathbf U\right)\mathbf s
=
\mathbf D_g\mathbf a^{\mathrm{ext}}.
}
\]

Defina

\[
\mathbf A=\mathbf I-\mathbf D_g\mathbf U,
\qquad
\mathbf b=\mathbf D_g\mathbf a^{\mathrm{ext}},
\]

e resolva inicialmente com `numpy.linalg.solve`.

Calcule e retorne:

\[
\mathcal R
=
\frac{\|\mathbf A\mathbf s-\mathbf b\|_2}
{\max(\|\mathbf b\|_2,\epsilon_{\mathrm{mach}})}
\]

e o número de condição

\[
\kappa_2(\mathbf A).
\]

Não forme explicitamente \(\mathbf A^{-1}\).

---

## 9. API mínima sugerida

É permitido ajustar nomes para manter a arquitetura limpa, mas a separação de
responsabilidades deve ser preservada.

### `src/acoustic_ms/multipoles.py`

```python
mode_count(lmax: int) -> int
mode_index(ell: int, m: int) -> int
mode_from_index(index: int) -> tuple[int, int]
modes(lmax: int) -> tuple[tuple[int, int], ...]
```

### `src/acoustic_ms/special.py`

```python
spherical_hankel1(ell: int, x: object, derivative: bool = False) -> object
spherical_harmonic(ell: int, m: int, theta: object, phi: object) -> object
cartesian_to_spherical(vector_xyz: object) -> tuple[float, float, float]
```

O retorno de `cartesian_to_spherical` deve ser `(r, theta, phi)`.

### `src/acoustic_ms/gaunt.py`

```python
gaunt_coefficient(
    n_source: int,
    m_source: int,
    q: int,
    mu: int,
    n_target: int,
    m_target: int,
) -> float
```

### `src/acoustic_ms/translation.py`

```python
separation_coefficient(
    n_target: int,
    m_target: int,
    n_source: int,
    m_source: int,
    k: float,
    target_position_xyz: object,
    source_position_xyz: object,
) -> complex

translation_matrix(
    k: float,
    target_position_xyz: object,
    source_position_xyz: object,
    target_lmax: int,
    source_lmax: int | None = None,
) -> numpy.ndarray
```

`translation_matrix` deve retornar uma matriz com shape

```text
(mode_count(target_lmax), mode_count(source_lmax))
```

e usar `source_lmax = target_lmax` por padrão.

### `src/acoustic_ms/scattering.py`

```python
rayleigh_scattering_coefficients(
    ka: float,
    f0: float,
    f1: float,
) -> numpy.ndarray
```

Retorna `[s0, s1]`; a repetição em \(m\) pertence à montagem da matriz
diagonal.

### `src/acoustic_ms/incident.py`

```python
nodal_standing_wave_coefficients(lmax: int) -> numpy.ndarray
```

### `src/acoustic_ms/solver.py`

Use um `dataclass` imutável ou estrutura equivalente contendo, no mínimo:

```python
coefficients          # complex ndarray, shape (n_particles, 4)
system_matrix         # complex ndarray, shape (4N, 4N)
right_hand_side       # complex ndarray, shape (4N,)
residual_relative     # float
condition_number      # float
modes                 # ordering used
```

API pública mínima:

```python
solve_rayleigh_nodal(
    positions_xyz: object,
    k: float,
    radius: float,
    f0: float,
    f1: float,
    lmax: int = 1,
) -> RayleighNodalSolution
```

Exporte apenas as funções públicas úteis em `src/acoustic_ms/__init__.py`.
Funções auxiliares internas não precisam compor a API pública.

---

## 10. Validação de entradas

Reutilize ou centralize validadores existentes quando isso reduzir duplicação,
sem alterar as APIs aprovadas.

O solver deve rejeitar:

- `k <= 0`;
- `radius <= 0`;
- valores não reais ou não finitos;
- `ka > 0.1`;
- `lmax != 1`;
- array de posições que não tenha shape `(N, 3)`;
- \(N<1\);
- centros fora do plano \(z=0\);
- centros coincidentes;
- partículas sobrepostas, \(r_{ij}<2a\);
- \(f_1\notin[-2,1]\).

Para a verificação do plano, use uma tolerância explícita, documentada e
proporcional à escala geométrica. Não use igualdade binária sem tolerância.

O operador de translação de baixo nível deve rejeitar centros coincidentes,
mas não precisa exigir que os centros estejam no plano nodal, pois seu teste
de reexpansão será tridimensional.

---

## 11. Testes obrigatórios

Crie, no mínimo:

```text
tests/test_multipoles.py
tests/test_special.py
tests/test_gaunt.py
tests/test_translation.py
tests/test_scattering.py
tests/test_incident.py
tests/test_solver.py
```

### 11.1 Indexação

- `mode_count(L) == (L + 1)**2`;
- ida e volta índice \(\leftrightarrow(\ell,m)\);
- ordem exata dos quatro modos para \(L_{\max}=1\);
- rejeição de índices inválidos.

### 11.2 Funções especiais

- \(h_\ell^{(1)}=j_\ell+iy_\ell\);
- derivadas conferidas por diferença finita ou identidade de recorrência;
- wronskiano

  \[
  j_\ell y_\ell'-j_\ell' y_\ell=\frac1{x^2};
  \]

- identidade \(Y_\ell^{-m}=(-1)^mY_\ell^{m*}\);
- ortonormalidade de alguns harmônicos por quadratura;
- ângulos cartesianos de vetores nos eixos e em um ponto geral.

### 11.3 Gaunt

Verifique valores conhecidos:

\[
\mathcal G(1,0;0,0;1,0)=\frac1{\sqrt{4\pi}},
\]

\[
\mathcal G(1,0;2,0;1,0)=\frac{\sqrt5}{5\sqrt\pi}.
\]

Além disso:

- teste zeros impostos pelas regras de seleção;
- compare alguns coeficientes não triviais com integração angular numérica
  independente.

### 11.4 Translação

- teste analítico

  \[
  \mathcal S_{10,10}=h_0^{(1)}(kd)+h_2^{(1)}(kd);
  \]

- teste direto do teorema de reexpansão descrito na Seção 5.3;
- convergência do erro quando \(L_{\mathrm{teste}}\) aumenta;
- shape e ordenação da matriz;
- rejeição de centros coincidentes.

### 11.5 Uma partícula

Para \(N_p=1\), apenas o coeficiente \((1,0)\) deve ser não nulo:

\[
\boxed{
s_{10}^{\mathrm{single}}
=s_1a_{10}
=if_1\sqrt{\frac{\pi}{3}}(ka)^3.
}
\]

Exija erro relativo menor que \(10^{-13}\) e resíduo compatível com precisão
de máquina.

### 11.6 Duas partículas

Considere partículas idênticas em

\[
\mathbf r_1=(-d/2,0,0),
\qquad
\mathbf r_2=(d/2,0,0).
\]

No truncamento \(L_{\max}=1\), a simetria deixa apenas \(s_{10}\), igual nas
duas partículas. Definindo

\[
C(kd)=h_0^{(1)}(kd)+h_2^{(1)}(kd),
\]

o resultado analítico é

\[
\boxed{
s_{10}^{\mathrm{pair}}
=
\frac{s_1a_{10}}
{1-s_1C(kd)}.
}
\]

Compare o solver global com essa expressão em pelo menos três casos,
incluindo:

- grande separação;
- separação intermediária;
- contato \(d=2a\).

Exija erro relativo menor que \(10^{-12}\).

No limite \(kd\ll1\), verifique também

\[
s_1C(kd)\longrightarrow
\frac{f_1}{2}\left(\frac ad\right)^3
\]

e, consequentemente,

\[
\frac{s_{10}^{\mathrm{pair}}}{s_{10}^{\mathrm{single}}}
\longrightarrow
\frac{1}
{1-\dfrac{f_1}{2}(a/d)^3}.
\]

### 11.7 Reespalhamento de todas as ordens

Para um caso com \(|s_1C|<1\), compare a solução direta com somas parciais da
série de Neumann:

\[
\mathbf s^{(P)}
=
\sum_{p=0}^{P}
(\mathbf D_g\mathbf U)^p
\mathbf D_g\mathbf a^{\mathrm{ext}}.
\]

Mostre no teste que o erro diminui com \(P\) e converge para a solução de
`numpy.linalg.solve`. Esse teste deixa explícito que \(L_{\max}=1\) não
significa apenas um ou dois eventos de espalhamento.

### 11.8 Três partículas

Use um triângulo escaleno no plano \(xy\) apenas como teste estrutural:

- resíduo relativo menor que \(10^{-12}\);
- covariância por permutação das partículas;
- invariância dos coeficientes \(s_{10}^{(i)}\) sob rotação rígida do cluster
  no plano;
- todos os modos proibidos pela simetria nodal devem permanecer nulos dentro
  da tolerância numérica.

Não calcule força e não apresente esse teste como resultado científico de três
corpos.

### 11.9 Regressão

Todos os 32 testes das T01 e T02 devem continuar passando, inclusive com
warnings promovidos a erro.

---

## 12. Dependências

Adicione às dependências de execução:

```toml
scipy
sympy
```

Mantenha NumPy. Não adicione bibliotecas de álgebra linear esparsa, Numba,
JAX, CuPy ou ferramentas de paralelização nesta tarefa.

O solver da T03 deve usar álgebra densa e CPU. Otimização pertence a uma etapa
posterior, depois da validação científica.

---

## 13. Script reprodutível de validação

Crie:

```text
scripts/validate_t03_solver.py
```

O script deve:

1. executar um caso de uma partícula;
2. executar três casos de duas partículas;
3. executar um triângulo escaleno;
4. registrar número de partículas, \(ka\), \(d/a\) quando aplicável,
   \(f_0\), \(f_1\), resíduo, número de condição e erro do benchmark quando
   houver solução analítica;
5. escrever

   ```text
   results/data/t03_solver_validation.csv
   ```

6. imprimir um resumo curto e determinístico.

Não gere gráfico nesta tarefa. O resultado principal é a validação do núcleo
matemático.

---

## 14. Documentação e rastreabilidade

Atualize:

- `README.md`;
- `TASKS.md`;
- `docs/CONVENTIONS.md`;
- `docs/DECISIONS.md`;
- `docs/HANDOFF.md`;
- o docstring de `src/acoustic_ms/__init__.py`.

Registre explicitamente:

- \(L_{\max}=1\) limita a ordem multipolar, não a quantidade de
  reespalhamentos;
- o solver contém \(4N_p\) incógnitas na base completa;
- a orientação do operador é alvo \(\leftarrow\) fonte;
- a matriz usa linhas-alvo e colunas-fonte;
- a solução usa \(\mathbf A=\mathbf I-\mathbf D_g\mathbf U\);
- a T03 resolve coeficientes de campo, não força;
- o solver de produção ainda não foi estendido para \(L_{\max}=3,5\).

Corrija também os dois detalhes cosméticos já identificados na auditoria:

- o título de `docs/HANDOFF.md` ainda começa como “T01 handoff”;
- o docstring de `src/acoustic_ms/__init__.py` ainda descreve somente a T01.

Não reescreva o histórico das T01 e T02; apenas organize títulos e acrescente a
seção da T03.

---

## 15. Critérios objetivos de aprovação

A T03 somente será aprovada se:

1. todos os testes antigos e novos passarem;
2. `pytest -q -W error` terminar sem falhas ou warnings;
3. o teste físico de reexpansão apresentar erro relativo \(<10^{-9}\);
4. o benchmark de uma partícula apresentar erro \(<10^{-13}\);
5. o benchmark analítico de duas partículas apresentar erro \(<10^{-12}\);
6. todos os casos do script apresentarem resíduo relativo \(<10^{-12}\);
7. o CSV puder ser regenerado sem alteração;
8. a orientação do operador estiver documentada e protegida por teste
   tridimensional fora dos eixos;
9. nenhuma força multibody tiver sido implementada ou alegada;
10. nenhuma API ou resultado aprovado das T01 e T02 tiver sido alterado.

---

## 16. Comandos finais obrigatórios

Execute, em ambiente limpo ou no ambiente virtual do projeto:

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python scripts/validate_t03_solver.py
python -m pytest -q -W error
git status --short
git diff --check
git diff --stat
```

Se o ambiente já possuir `.venv`, use o interpretador correspondente.

---

## 17. Relatório final para auditoria

Ao terminar, responda com:

1. resumo do que foi implementado;
2. lista dos arquivos criados e alterados;
3. equações efetivamente implementadas;
4. convenção de deslocamento usada na translação;
5. resultado do teste direto de reexpansão;
6. resultado dos benchmarks \(N_p=1\) e \(N_p=2\);
7. resíduos e números de condição;
8. resultado completo do `pytest -q -W error`;
9. resumo do `git diff --stat`;
10. limitações e dúvidas científicas;
11. confirmação explícita de que nenhuma força multibody foi implementada.

Inclua em `docs/HANDOFF.md` as mesmas informações essenciais.

Não esconda tolerâncias relaxadas, warnings, testes pulados, condicionamento
ruim ou diferenças em relação a esta especificação.
