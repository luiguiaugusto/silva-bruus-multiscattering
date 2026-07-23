# T02 — Fórmula corrigida de duas partículas e reprodução da Figura 2

## Instrução de uso

Copie esta especificação integralmente para o Codex CLI aberto na raiz do
repositório `silva-bruus-multiscattering`.

Esta tarefa parte do estado aprovado da T01. Antes de alterar qualquer arquivo,
leia `AGENTS.md`, `README.md`, `TASKS.md`, `docs/CONVENTIONS.md`,
`docs/DECISIONS.md`, `docs/HANDOFF.md`, o código em `src/acoustic_ms/` e os
testes existentes.

Não faça `commit` nem `push` nesta tarefa, a menos que o usuário solicite isso
separadamente.

---

## 1. Objetivo

Implemente a expressão analítica corrigida da força de interação acústica entre
duas esferas idênticas no plano nodal, apresentada nas Eqs. (30a)-(30d) do
artigo de 2026.

Em seguida:

1. compare essa expressão com a força de Silva-Bruus já implementada na T01;
2. reproduza numericamente a Figura 2 do artigo;
3. gere os dados e a figura de maneira reprodutível;
4. crie testes que recuperem os quatro erros máximos publicados.

A T02 produz o **Modelo B de duas partículas**. Ela ainda é um benchmark
conhecido e não constitui o modelo de múltiplo espalhamento para \(N\geq3\).

---

## 2. Escopo físico

Preserve todas as convenções da T01:

\[
p(\mathbf r,t)=\operatorname{Re}
\{p(\mathbf r)e^{-i\omega t}\},
\qquad
E_0=\frac{\rho_0|v_0|^2}{4}.
\]

Considere:

- fluido ideal, homogêneo, invíscido e ilimitado;
- acústica linear e harmônica;
- duas esferas compressíveis, idênticas e não sobrepostas;
- partículas contidas no plano nodal de pressão;
- regime de Rayleigh;
- posições fixas;
- força dirigida ao longo da linha de centros;
- ausência de paredes, viscosidade, streaming e interação hidrodinâmica.

Para a força sobre a partícula-prova \(i\) devida à partícula-fonte \(j\),
continue usando

\[
\mathbf d_{ij}=\mathbf r_i-\mathbf r_j,
\qquad
d=\|\mathbf d_{ij}\|,
\qquad
\widehat{\mathbf d}_{ij}=\frac{\mathbf d_{ij}}{d}.
\]

O vetor \(\widehat{\mathbf d}_{ij}\) aponta da fonte para a prova. Portanto:

- componente radial positiva: repulsão;
- componente radial negativa: atração.

A condição de não sobreposição permanece

\[
d\geq 2a.
\]

Não trate \(ka\), \(kd\) e \(a/d\) como parâmetros independentes. Sempre vale

\[
kd=(ka)\frac{d}{a},
\qquad
\chi=\frac{a}{d}=\frac{ka}{kd}.
\]

---

## 3. Equações obrigatórias

Defina

\[
x=kd,
\qquad
\chi=\frac{a}{d}.
\]

A força corrigida é

\[
\boxed{
\mathbf F_{i\leftarrow j}^{\mathrm{corr}}
=
4\pi a^2 f_1^2E_0\chi^4
\left[
\frac{3A_0}{D_0}
\bigl(\cos x+x\sin x\bigr)
+
\frac{A_2}{5D_0}x^2\cos x
\right]
\widehat{\mathbf d}_{ij}
}
\]

com

\[
A_0=
4(77-25f_1+2f_1^2)
+90f_1(-11+2f_1)\chi^7
+5250(-7+f_1)f_1\chi^{11}
+33075f_1^2\chi^{18},
\]

\[
A_2=
-20(77-25f_1+2f_1^2)
-162f_1(-11+2f_1)\chi^7
-24250(-7+f_1)f_1\chi^{11}
+20025f_1^2\chi^{18},
\]

e

\[
\begin{aligned}
D_0={}&-616
+4f_1\left(50+77\chi^3+1485\chi^7+22050\chi^{11}\right)
\\
&+f_1^3\chi^3
\left(8+324\chi^7+5850\chi^{11}+73575\chi^{18}\right)
\\
&-2f_1^2
\left[
8+50\chi^3
+9\chi^7
\left(
60+\chi^3
\left(
99+700\chi+2275\chi^4+17850\chi^8
\right)
\right)
\right].
\end{aligned}
\]

Implemente os parênteses exatamente como escritos. Não simplifique
manualmente os polinômios no código se isso tornar a expressão mais difícil de
auditar contra o artigo.

A força de referência da T01 é

\[
\boxed{
\mathbf F_{i\leftarrow j}^{\mathrm{SB}}
=
4\pi a^2E_0f_1^2\chi^4
\left[
-\frac32\bigl(\cos x+x\sin x\bigr)
+\frac12x^2\cos x
\right]
\widehat{\mathbf d}_{ij}
}.
\]

Use obrigatoriamente a função já existente na T01 para calcular
\(\mathbf F^{\mathrm{SB}}\). Não crie uma segunda implementação da fórmula de
Silva-Bruus.

---

## 4. Limites analíticos que devem ser preservados

### 4.1 Contraste de densidade fraco

Quando \(f_1\to0\),

\[
\frac{A_0}{D_0}\to-\frac12,
\qquad
\frac{A_2}{D_0}\to\frac52,
\]

e a expressão corrigida tende à expressão de Silva-Bruus.

Não valide esse limite comparando apenas as forças em \(f_1=0\), pois ambas
seriam trivialmente nulas devido ao prefator \(f_1^2\). Teste diretamente as
razões \(A_0/D_0\) e \(A_2/D_0\), além de comparar as forças para um
\(f_1\) pequeno, mas não nulo.

### 4.2 Grande separação geométrica

Quando \(\chi=a/d\to0\), para qualquer \(f_1\) admissível,

\[
\frac{A_0}{D_0}\to-\frac12,
\qquad
\frac{A_2}{D_0}\to\frac52,
\]

e novamente

\[
\mathbf F^{\mathrm{corr}}\to\mathbf F^{\mathrm{SB}}.
\]

### 4.3 Pequeno \(kd\)

Para \(kd\ll1\),

\[
\mathbf F^{\mathrm{corr}}
=
12\pi a^2f_1^2E_0
\left(\frac ad\right)^4
\frac{A_0}{D_0}\widehat{\mathbf d}
+O\!\left((kd)^2\right).
\]

No regime estudado no artigo, \(A_0/D_0<0\), logo a força próxima ao contato é
atrativa.

---

## 5. API mínima

Crie

```text
src/acoustic_ms/corrected_pair.py
```

com uma API pública mínima equivalente a:

```python
corrected_pair_coefficients(
    f1: float,
    separation_ratio: float,
) -> tuple[float, float, float]
```

Retorna `(A0, A2, D0)`.

```python
corrected_nodal_pair_force_magnitude(
    k: float,
    radius: float,
    distance: float,
    energy_density: float,
    f1: float,
) -> float
```

Apesar de o nome seguir a API existente da T01, o retorno é uma **componente
radial com sinal**, em newtons, e não um módulo necessariamente positivo. Isso
deve estar explícito na docstring.

```python
corrected_nodal_pair_force_on_probe(
    probe_xy: object,
    source_xy: object,
    k: float,
    radius: float,
    energy_density: float,
    f1: float,
) -> numpy.ndarray
```

Retorna a força vetorial 2D sobre a prova.

```python
corrected_nodal_pair_forces(
    position_1_xy: object,
    position_2_xy: object,
    k: float,
    radius: float,
    energy_density: float,
    f1: float,
) -> tuple[numpy.ndarray, numpy.ndarray]
```

Retorna as forças iguais e opostas sobre as duas partículas.

Exporte as novas funções em `src/acoustic_ms/__init__.py`.

Reutilize as regras de validação da T01 e evite duplicação desnecessária. Uma
refatoração pequena das rotinas internas de validação é permitida, desde que:

- a API pública da T01 permaneça compatível;
- todos os 16 testes existentes continuem passando;
- o comportamento físico da T01 não seja alterado.

Para `corrected_pair_coefficients`, aceite somente escalares reais e finitos,
\(-2\leq f_1\leq1\) e

\[
0<\chi\leq\frac12.
\]

Na API de força, obtenha \(\chi\) de `radius / distance`; não o receba como
parâmetro independente.

---

## 6. Definição do erro da Figura 2

Para reproduzir os números publicados, use a solução corrigida como
referência:

\[
\boxed{
\varepsilon_{\mathrm{SB}}(\%)
=
100
\frac{
\left|F^{\mathrm{corr}}-F^{\mathrm{SB}}\right|
}{
\left|F^{\mathrm{corr}}\right|
}
}.
\]

Não use \(|F^{\mathrm{SB}}|\) no denominador e não use o erro relativo
simétrico nesta figura.

Como as duas forças são colineares no problema de duas partículas, a expressão
pode ser calculada com suas componentes radiais com sinal. O erro é indefinido
quando \(F^{\mathrm{corr}}=0\); não esconda essa situação adicionando um
epsilon arbitrário ao denominador.

---

## 7. Reprodução da Figura 2

Crie

```text
scripts/reproduce_figure_2.py
```

O script deve:

1. usar `ka = 0.1`;
2. percorrer \(0.2\leq kd\leq0.3\) com resolução suficiente para curvas
   suaves, por exemplo 501 pontos;
3. usar \(f_1=0.1,\ 0.4,\ 0.8,\ 1.0\);
4. calcular, para cada ponto,
   \[
   \chi=\frac{ka}{kd};
   \]
5. chamar as funções científicas do pacote, sem recodificar as equações no
   script;
6. calcular o erro percentual com
   \(F^{\mathrm{corr}}\) no denominador;
7. salvar os dados em
   `results/data/figure_2_relative_error.csv`;
8. salvar a figura em
   `results/figures/figure_2_relative_error.png`.

O CSV deve conter pelo menos:

```text
ka,kd,separation_ratio,f1,corrected_force,silva_bruus_force,relative_error_percent
```

Use uma escala SI simples e documentada para avaliar as forças, por exemplo
`radius = 1`, `k = ka / radius` e `energy_density = 1`. O erro é independente
dessa escolha porque o prefator comum cancela.

A figura deve conter:

- eixo horizontal \(kd\);
- eixo vertical `Relative error (%)`;
- uma curva para cada \(f_1\);
- legenda;
- limites horizontais de 0.20 a 0.30;
- resolução adequada para leitura e uso no relatório.

O artigo escreve o intervalo como \(2ka<kd<0.3\), mas os erros máximos
publicados correspondem ao limite de contato

\[
kd\to2ka=0.2,
\qquad
\chi\to\frac12.
\]

O código atual admite \(d=2a\). Portanto, use \(kd=0.2\) como ponto-limite
para o teste de regressão e para recuperar os máximos, deixando isso explícito
na documentação e na legenda ou descrição do resultado. Não interprete esse
ponto como uma separação com folga positiva entre as superfícies.

Adicione o Matplotlib como dependência opcional claramente documentada, sem
obrigar o núcleo de cálculo a importar Matplotlib.

---

## 8. Valores independentes de referência

Use os valores abaixo como testes de regressão. Eles foram calculados
independentemente a partir das Eqs. (30a)-(30d), com

\[
ka=0.1,\qquad kd=0.2,\qquad \chi=0.5.
\]

### 8.1 Razões dos coeficientes

| \(f_1\) | \(A_0/D_0\) | \(A_2/D_0\) |
|---:|---:|---:|
| 0.1 | -0.506355807146221 | 2.537084546937989 |
| 0.4 | -0.527267939835848 | 2.659985720443811 |
| 0.8 | -0.560671850474030 | 2.858882625216557 |
| 1.0 | -0.580568297988085 | 2.978697291068864 |

### 8.2 Erros no limite de contato

| \(f_1\) | erro calculado (%) | valor arredondado do artigo (%) |
|---:|---:|---:|
| 0.1 | 1.252519728707 | 1.25 |
| 0.4 | 5.160511340274 | 5.16 |
| 0.8 | 10.798343941865 | 10.80 |
| 1.0 | 13.848266387733 | 13.84 aproximadamente |

O último valor pode aparecer como 13.85% ao arredondar numericamente a duas
casas. O critério científico é concordar com a curva e com a afirmação
“aproximadamente 13.84%”, e não forçar arredondamento artificial no código.

Para um teste adicional, em \(kd=0.3\), \(ka=0.1\) e \(f_1=1\), o erro deve
ser aproximadamente

\[
2.206083641258\%.
\]

---

## 9. Testes obrigatórios

Crie

```text
tests/test_corrected_pair.py
```

e cubra, no mínimo:

1. valores independentes de \(A_0/D_0\) e \(A_2/D_0\) da tabela;
2. os quatro erros no limite \(kd=0.2\);
3. o ponto adicional \(kd=0.3,\ f_1=1\);
4. \(A_0/D_0=-1/2\) e \(A_2/D_0=5/2\) para \(f_1=0\);
5. convergência da força corrigida para Silva-Bruus quando
   \(f_1\to0\), usando contraste pequeno não nulo;
6. convergência da força corrigida para Silva-Bruus quando
   \(a/d\to0\);
7. limite assintótico de pequeno \(kd\);
8. força atrativa no regime próximo usado na Figura 2;
9. ação e reação;
10. invariância por translação;
11. covariância por rotação no plano;
12. força nula para `energy_density = 0`;
13. força nula para `f1 = 0`;
14. rejeição de sobreposição \(d<2a\);
15. rejeição de posições coincidentes;
16. rejeição de escalares não finitos e de valores fora dos domínios físicos;
17. diminuição monótona do erro com \(kd\) entre 0.2 e 0.3 para cada um dos
    quatro \(f_1\) da Figura 2;
18. preservação integral dos testes da T01.

Os testes dos números publicados devem chamar as funções públicas do pacote e
compará-las com constantes de referência. Não escreva, dentro do teste, uma
cópia algébrica da mesma fórmula usada na implementação.

Use tolerâncias numéricas justificadas. Para os valores completos fornecidos
acima, uma tolerância relativa da ordem de \(10^{-11}\) a \(10^{-12}\) é
adequada em precisão dupla. Para os valores arredondados do artigo, aceite a
incerteza do arredondamento decimal.

---

## 10. Arquivos a criar ou atualizar

Crie:

- `src/acoustic_ms/corrected_pair.py`;
- `tests/test_corrected_pair.py`;
- `scripts/reproduce_figure_2.py`;
- `results/data/figure_2_relative_error.csv`;
- `results/figures/figure_2_relative_error.png`.

Atualize somente quando necessário:

- `src/acoustic_ms/__init__.py`;
- `pyproject.toml`;
- `README.md`;
- `TASKS.md`;
- `docs/CONVENTIONS.md`;
- `docs/DECISIONS.md`;
- `docs/HANDOFF.md`;
- `.gitignore`, apenas se algum artefato transitório novo precisar ser
  ignorado.

Não altere os PDFs de referência.

---

## 11. Fora do escopo

Não implemente nesta tarefa:

- o sistema multipolar que originou a fórmula fechada;
- harmônicos esféricos;
- funções de Bessel, Neumann ou Hankel;
- matrizes T;
- operadores de translação;
- solução de sistema linear acoplado;
- clusters com \(N\geq3\);
- soma de vários pares corrigidos;
- série de Neumann;
- Figura 3, dados experimentais ou comparação com poliestireno;
- varreduras gerais em \(ka\), \(f_1\) ou geometrias;
- otimizações, Numba, GPU ou paralelização.

Em particular, não alegue que a T02 implementa uma “solução completa de
múltiplo espalhamento”. Ela implementa somente a **expressão analítica
truncada em quinta ordem publicada para duas partículas**.

---

## 12. Verificação e critérios de conclusão

Execute, em ambiente limpo ou no ambiente virtual existente:

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python scripts/reproduce_figure_2.py
python -m pytest -q -W error
```

Adapte apenas o nome do extra opcional caso escolha outro nome claro no
`pyproject.toml`.

A tarefa estará concluída somente se:

1. todos os 16 testes antigos continuarem passando;
2. todos os testes novos passarem;
3. a fórmula usar exatamente as Eqs. (30a)-(30d);
4. a API vetorial preservar a convenção da T01;
5. os quatro erros máximos forem recuperados;
6. o erro usar \(F^{\mathrm{corr}}\) no denominador;
7. o CSV for gerado pelas funções do pacote;
8. a figura for gerada a partir do CSV ou dos mesmos dados calculados, sem
   valores hard-coded;
9. `docs/HANDOFF.md` registrar equações, decisões, comandos, resultados e
   limitações;
10. `TASKS.md` marcar somente a T02 como concluída;
11. nenhuma implementação da T03 ou de \(N\geq3\) for adicionada.

---

## 13. Relatório obrigatório ao final

Ao terminar, mostre:

```bash
git status --short
git diff --stat
python -m pytest -q
```

e apresente um resumo com:

```markdown
# Relatório T02

## Arquivos criados ou alterados

## Equações implementadas

## Definição de erro utilizada

## Testes executados

## Valores reproduzidos da Figura 2

## Caminhos do CSV e da figura

## Limitações ou dúvidas

## Resumo do diff
```

Se surgir qualquer conflito entre esta especificação e um arquivo antigo do
repositório, preserve as convenções aprovadas da T01 e relate o conflito. Não
invente uma nova convenção física.

---

## 14. Referência primária

Use como fonte das equações o PDF do artigo de 2026 já presente no repositório:

```text
papers/Acoustic_Interaction_Force (1)(2).pdf
```

Pontos de conferência:

- Eqs. (30a)-(30d): força corrigida e coeficientes;
- Eqs. (31a)-(31b): limite \(f_1\to0\);
- Eq. (33): limite de pequeno \(kd\);
- Eq. (38): força de Silva-Bruus;
- Figura 2: benchmark de erro relativo.

O código existente da T01 é a referência executável para a Eq. (38).
