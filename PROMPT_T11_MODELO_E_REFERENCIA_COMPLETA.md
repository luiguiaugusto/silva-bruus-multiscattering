# Prompt para execução da T11 - Referência multipolar completa, Modelo E

Você está trabalhando no repositório:

```text
silva-bruus-multiscattering
```

Execute a **T11 - construção e validação do Modelo E**, uma referência
multipolar completa para esferas fluidas idênticas em um plano nodal.

Esta tarefa começa após o fechamento da T10. A base esperada da `main` é:

```text
153403c feat: add exact Mie scattering coefficients
```

O objetivo científico é combinar, pela primeira vez no projeto:

1. os coeficientes exatos de Mie implementados na T10;
2. um sistema global de múltiplo espalhamento;
3. a força de radiação multipolar completa sobre cada esfera;
4. a separação entre as contribuições `external-scattered` e
   `scattered-scattered`.

O Modelo E será o oráculo físico usado posteriormente na T12 para avaliar os
Modelos A e D. A T11 deve construir e validar esse oráculo, mas **não deve
executar ainda a campanha científica de 20-30 casos da T12**.

---

## 1. Auditoria inicial obrigatória

Antes de editar:

1. Leia integralmente:

   - `AGENTS.md`;
   - `README.md`;
   - `TASKS.md`;
   - `docs/CONVENTIONS.md`;
   - `docs/DECISIONS.md`;
   - `docs/HANDOFF.md`;
   - `TAREFA_T07_MODELO_D_CONVERGENCIA_MULTIPOLAR.md`;
   - `TAREFA_T09_FUNDAMENTACAO_ANALITICA_RHO1.md`;
   - `TAREFA_T10_COEFICIENTES_EXATOS_MIE.md`;
   - `src/acoustic_ms/incident.py`;
   - `src/acoustic_ms/multipoles.py`;
   - `src/acoustic_ms/translation.py`;
   - `src/acoustic_ms/multipolar_scattering.py`;
   - `src/acoustic_ms/multipolar_solver.py`;
   - `src/acoustic_ms/model_d.py`;
   - `src/acoustic_ms/force.py`;
   - `src/acoustic_ms/mie_scattering.py`;
   - os testes associados a esses módulos.

2. Consulte as referências:

   - `papers/Acoustic_Interaction_Force (1)(2).pdf`, especialmente as
     Eqs. (11), (16), (21) e (22);
   - `papers/2014-silva-bruus(2).pdf`, especialmente a aproximação pairwise;
   - Lopes, Azarpeyvand e Silva, *IEEE Transactions on Ultrasonics,
     Ferroelectrics, and Frequency Control* **63**, 186-197 (2016),
     DOI `10.1109/TUFFC.2015.2494693`, especialmente as Eqs. (20), (34),
     (35), (43) e (44).

   As fórmulas necessárias da terceira referência estão explicitadas neste
   prompt. Não dependa de acesso à internet para implementá-las.

3. Execute:

   ```bash
   git status --short
   git rev-parse --short HEAD
   git log -1 --oneline
   ```

4. O diretório de trabalho deve estar limpo e o `HEAD` deve ser `153403c`.
   Se estiver apenas desatualizado e sem mudanças locais, execute:

   ```bash
   git pull --ff-only origin main
   ```

5. Se houver mudanças locais, conflitos ou uma base diferente após o `pull`,
   não descarte nem sobrescreva nada. Pare e informe o estado encontrado.

6. Registre os hashes dos artefatos já versionados sob `results/data/` e
   `results/figures/`. Eles devem permanecer byte a byte inalterados.

7. Instale o projeto, se necessário, e execute a suíte-base:

   ```bash
   .venv/bin/python -m pip install -e ".[dev,plot]"
   .venv/bin/python -m pytest -q -W error
   ```

O resultado-base esperado é:

```text
258 passed
```

Se a suíte-base falhar, não comece a implementação antes de diagnosticar e
relatar a falha.

---

## 2. Definição científica do Modelo E

O Modelo E deve permanecer especializado ao domínio físico atual:

- fluido hospedeiro ideal, invíscido e ilimitado;
- esferas fluidas, compressíveis, homogêneas, idênticas e sem perdas;
- dependência temporal \(e^{-i\omega t}\);
- onda estacionária nodal \(\widetilde p_{\mathrm{ext}}=\sin(kz)\);
- centros no plano \(z=0\);
- partículas não sobrepostas;
- T-matrix diagonal exata da T10;
- múltiplo espalhamento global resolvido em um truncamento multipolar comum;
- força de radiação completa calculada com os coeficientes incidentes efetivos.

Não renomeie nem altere os Modelos A-D. O Modelo E deve ser uma API nova.

### 2.1 Coeficientes externos, incidentes efetivos e espalhados

Para cada partícula \(i\), denote por

\[
a_{nm}^{(i)}
\]

os coeficientes da onda externa, por

\[
b_{nm}^{(i)}
\]

os coeficientes do campo incidente efetivo e por

\[
d_{nm}^{(i)}=s_n^{\mathrm{Mie}}b_{nm}^{(i)}
\]

os coeficientes do campo espalhado pela própria partícula.

Use a letra \(d_{nm}\), ou outro nome não ambíguo no código, para não
confundir o vetor de coeficientes espalhados com o coeficiente escalar
\(s_n^{\mathrm{Mie}}\).

O campo incidente efetivo é

\[
b_i=a_i+\sum_{j\ne i}U_{ij}d_j,
\qquad
d_i=D b_i,
\]

onde \(D\) é a T-matrix diagonal formada pelos coeficientes exatos da T10 e
\(U_{ij}\) é a translação `target <- source` já adotada no projeto.

Portanto, o sistema principal do Modelo E deve ser montado na forma

\[
\boxed{
(\mathbf I-\mathbf U\mathbf D)\mathbf b=\mathbf a
}
\]

e resolvido com `numpy.linalg.solve`, nunca por inversão explícita.

Após a solução:

\[
\boxed{
\mathbf d=\mathbf D\mathbf b.
}
\]

Essa formulação em termos do campo incidente efetivo é deliberadamente
diferente da API de saída do Modelo D e coincide com a formulação de Lopes,
Azarpeyvand e Silva.

Como verificação independente do sistema, implemente em teste a forma
equivalente para os coeficientes espalhados:

\[
(\mathbf I-\mathbf D\mathbf U)\mathbf d=\mathbf D\mathbf a.
\]

As duas formas devem concordar, mas a segunda não deve ser o caminho de
produção do Modelo E.

### 2.2 Coeficientes exatos de Mie

Use exclusivamente:

```python
mie_scattering_coefficients_from_contrasts(...)
```

da T10 para formar \(D\). Não copie a fórmula de Mie para outro módulo.

Preserve:

- \(f_0<1\);
- \(-2<f_1\le1\);
- o limite rígido analítico em \(f_1=1\);
- o caso de correspondência material;
- a independência de \(m\) para cada ordem \(n\).

O Modelo E não pode chamar
`rayleigh_multipolar_scattering_coefficients` em seu caminho de produção.
Coeficientes de Rayleigh podem aparecer somente em testes de redução
assintótica.

### 2.3 Base completa e simetria planar

O ordenamento completo deve continuar sendo:

\[
\operatorname{index}(n,m)=n^2+n+m.
\]

Implemente:

- uma solução com todos os \((n,m)\) até \(L_{\max}\);
- a redução opcional pela simetria planar \(n+m\) ímpar;
- coeficientes completos de saída, com zeros exatos nos canais eliminados;
- uma comparação automática entre base completa e base reduzida.

Não elimine genericamente todas as ordens pares. A regra correta para um
cluster planar genérico no nó é:

\[
\boxed{n+m\ \text{ímpar}.}
\]

Não remova um modo apenas porque um coeficiente de Mie é numericamente
pequeno. Transparência exata pode ser tratada como um caso especial, mas não
use limiar empírico para podar a base.

### 2.4 Truncamento da força

A força acopla ordens adjacentes \(n\) e \(n+1\). Para uma solução disponível
até \(L_{\max}=M\), some:

\[
n=0,\ldots,M-1.
\]

Nunca invente, extrapole ou preencha silenciosamente coeficientes de ordem
\(M+1\).

No problema nodal, \(M=1\) não contém ordens adjacentes suficientes para uma
força transversal não trivial. A API principal do Modelo E deve exigir:

\[
\boxed{L_{\max}\ge2.}
\]

A convergência do Modelo E deve ser verificada variando \(L_{\max}\), e não
presumida a partir do valor usado no Modelo D.

---

## 3. Força multipolar completa

### 3.1 Funcional de força

Implemente uma rotina genérica que recebe:

- coeficientes incidentes \(b_{nm}\) relativos ao centro da partícula;
- coeficientes escalares \(s_n\) até \(M\);
- \(k\);
- a densidade de energia do projeto \(E_0\);
- o ordenamento completo de modos.

Defina:

\[
\Gamma_n
=
s_n+s_{n+1}^{*}+2s_ns_{n+1}^{*}.
\]

O termo \(2s_ns_{n+1}^{*}\) é o termo de recuo quadrático associado ao campo
espalhado pela própria esfera na fórmula exata de força. Ele deve ser
mantido.

Com a normalização de energia do projeto, calcule:

\[
\boxed{
F_x+iF_y
=
\frac{iE_0}{k^2}
\sum_{n=0}^{M-1}\sum_{m=-n}^{n}
\sqrt{
\frac{(n+m+1)(n+m+2)}
{(2n+1)(2n+3)}
}
\left[
\Gamma_n b_{nm}b_{n+1,m+1}^{*}
+
\Gamma_n^{*}b_{n,-m}^{*}b_{n+1,-m-1}
\right]
}
\]

e

\[
\boxed{
F_z
=
\frac{2E_0}{k^2}
\operatorname{Im}
\sum_{n=0}^{M-1}\sum_{m=-n}^{n}
\sqrt{
\frac{(n-m+1)(n+m+1)}
{(2n+1)(2n+3)}
}
\Gamma_n b_{nm}b_{n+1,m}^{*}.
}
\]

Recupere:

```python
Fx = np.real(Fx_plus_iFy)
Fy = np.imag(Fx_plus_iFy)
```

e retorne um vetor real \((F_x,F_y,F_z)\).

### 3.2 Conversão da normalização de energia

A referência de Lopes, Azarpeyvand e Silva usa

\[
E_{\mathrm{LAS}}=\frac{p_0^2}{2\rho_0c_0^2}.
\]

O projeto usa

\[
E_0=\frac{\rho_0|v_0|^2}{4}
=
\frac{p_0^2}{4\rho_0c_0^2}.
\]

Logo:

\[
E_{\mathrm{LAS}}=2E_0.
\]

Os fatores nas duas fórmulas acima já incorporam essa conversão. Não introduza
um fator ajustado empiricamente. A normalização deve ser validada contra o
limite Rayleigh do Modelo D e contra o oráculo de tensor de radiação descrito
adiante.

### 3.3 Força externa e força de interação

Defina o funcional quadrático de força da esfera por:

\[
\mathcal F[q]
=
\text{força calculada pelas fórmulas acima usando }b=q.
\]

Para cada partícula:

\[
c=b-a,
\]

onde \(c\) é o campo incidente sobre a partícula proveniente das demais
partículas.

Calcule:

\[
\mathbf F_{\mathrm{total}}=\mathcal F[b],
\]

\[
\mathbf F_{\mathrm{external}}=\mathcal F[a],
\]

\[
\boxed{
\mathbf F_{\mathrm{int}}
=
\mathbf F_{\mathrm{total}}-\mathbf F_{\mathrm{external}}.
}
\]

Essa definição corresponde à Eq. (44) de Lopes, Azarpeyvand e Silva.

### 3.4 Decomposição external-scattered/scattered-scattered

Como \(\mathcal F\) é um funcional quadrático:

\[
\boxed{
\mathbf F_{\mathrm{ss}}=\mathcal F[c]
}
\]

é a contribuição quadrática do campo incidente espalhado pelas demais
partículas, e

\[
\boxed{
\mathbf F_{\mathrm{ext-sc}}
=
\mathcal F[b]-\mathcal F[a]-\mathcal F[c]
}
\]

é a parte cruzada entre a onda externa e o campo incidente espalhado.

Exija numericamente:

\[
\boxed{
\mathbf F_{\mathrm{int}}
=
\mathbf F_{\mathrm{ext-sc}}
+
\mathbf F_{\mathrm{ss}}.
}
\]

Essa é a decomposição relevante para comparar o Modelo E com a força parcial
do Modelo D.

Não confunda:

- \(\mathbf F_{\mathrm{ss}}=\mathcal F[c]\), que é quadrática no campo
  incidente proveniente das outras partículas;
- \(2s_ns_{n+1}^{*}\), que faz parte de \(\Gamma_n\) e representa o recuo do
  campo espalhado pela própria esfera.

Os dois efeitos são distintos e ambos devem permanecer no Modelo E.

### 3.5 Interpretação no plano nodal

Para centros exatamente em \(z=0\):

- a força externa isolada deve ser nula no truncamento compatível com a
  simetria;
- a força de interação relevante para o projeto é o vetor \(xy\);
- \(F_z\) deve ser mantido e auditado como diagnóstico de simetria;
- não force \(F_z=0\) manualmente;
- não imponha soma global das forças igual a zero em geometrias genéricas.

---

## 4. APIs e arquivos científicos requeridos

Crie, preferencialmente:

```text
src/acoustic_ms/mie_multiparticle.py
src/acoustic_ms/complete_force.py
src/acoustic_ms/model_e.py
```

A divisão exata pode ser ajustada se houver uma razão arquitetural clara, mas
as três responsabilidades devem permanecer separadas:

1. solução global exata de Mie;
2. funcional genérico de força completa;
3. orquestração e decomposição do Modelo E.

### 4.1 Solução global

Crie uma dataclass equivalente a:

```python
@dataclass(frozen=True)
class MieMultiparticleSolution:
    effective_incident_coefficients: np.ndarray
    scattered_coefficients: np.ndarray
    external_coefficients: np.ndarray
    scattering_by_ell: np.ndarray
    system_matrix: np.ndarray
    right_hand_side: np.ndarray
    residual_relative: float
    condition_number: float
    modes: tuple[tuple[int, int], ...]
    active_modes: tuple[tuple[int, int], ...]
    active_mode_indices: tuple[int, ...]
    lmax: int
    used_planar_symmetry: bool
```

Os arrays de coeficientes por partícula devem usar sempre o ordenamento
completo. A matriz do sistema pode permanecer na base ativa.

Crie uma função pública semelhante a:

```python
solve_mie_multiparticle_nodal(
    positions_xyz,
    k,
    radius,
    f0,
    f1,
    lmax,
    *,
    use_planar_symmetry=True,
) -> MieMultiparticleSolution
```

Valide:

- `positions_xyz` com forma `(N, 3)`;
- \(N\ge1\);
- centros em \(z=0\), dentro da tolerância já documentada;
- \(k>0\);
- \(a>0\);
- não sobreposição;
- `lmax` inteiro e positivo;
- domínios de \(f_0,f_1\) da T10;
- finitude de todos os argumentos.

Não imponha \(ka\le0.1\) dentro do solver genérico de Mie. A campanha desta
tarefa continuará restrita a \(ka\le0.1\), mas a API exata não deve herdar
silenciosamente a restrição do solver de Rayleigh.

### 4.2 Funcional de força

Crie uma função pública de baixo nível semelhante a:

```python
complete_radiation_force_from_bsc(
    incident_coefficients,
    scattering_by_ell,
    k,
    energy_density,
) -> np.ndarray
```

Ela deve:

- aceitar um vetor completo de coeficientes de uma partícula;
- inferir ou validar \(L_{\max}\);
- usar apenas ordens adjacentes disponíveis;
- retornar exatamente três componentes reais;
- ser invariável sob uma fase global aplicada aos coeficientes;
- escalar linearmente com `energy_density`;
- retornar zero exato quando `energy_density == 0`;
- rejeitar entradas incompatíveis ou não finitas.

### 4.3 Resultado do Modelo E

Crie uma dataclass equivalente a:

```python
@dataclass(frozen=True)
class ModelENodalResult:
    solution: MieMultiparticleSolution
    total_forces_xyz: np.ndarray
    external_forces_xyz: np.ndarray
    interaction_forces_xyz: np.ndarray
    external_scattered_forces_xyz: np.ndarray
    scattered_scattered_forces_xyz: np.ndarray
    incoming_scattered_coefficients: np.ndarray
    decomposition_residual_relative: float
```

Crie:

```python
solve_model_e_nodal(
    positions_xyz,
    k,
    radius,
    energy_density,
    f0,
    f1,
    lmax,
    *,
    use_planar_symmetry=True,
) -> ModelENodalResult
```

Exponha propriedades ou auxiliares `forces_xy` apenas se não ocultarem que a
implementação calcula três componentes.

Atualize `src/acoustic_ms/__init__.py` apenas depois que as APIs estiverem
validadas.

---

## 5. Oráculo independente por tensor de radiação

A fórmula analítica de força não pode validar a si própria. Implemente um
oráculo independente, usado somente em testes ou validação, que integre
numericamente o tensor de radiação sobre uma superfície esférica de controle
centrada em uma partícula.

Ele não pode:

- chamar `complete_radiation_force_from_bsc`;
- reutilizar as somas com \(\Gamma_n\);
- ajustar um prefator para concordar com a rotina de produção;
- usar diferenças finitas da energia de interação do Modelo E.

### 5.1 Campo local

Reconstrua a velocidade potencial adimensional local:

\[
\psi(\mathbf r)
=
\sum_{n,m}
\left[
b_{nm}j_n(kr)
+
s_nb_{nm}h_n^{(1)}(kr)
\right]
Y_n^m(\theta,\phi).
\]

Defina:

\[
\mathbf g=\frac{1}{k}\nabla\psi.
\]

As derivadas podem ser analíticas em coordenadas esféricas ou avaliadas por
um esquema numérico de ordem alta com estudo explícito de refinamento. Prefira
derivadas analíticas. Use quadratura de Gauss-Legendre em
\(\mu=\cos\theta\) e uma malha periódica uniforme em \(\phi\), evitando os
polos.

### 5.2 Tensor e integração

Com a normalização energética do projeto:

\[
\frac{\overline{\mathbf S}}{E_0}
=
-
\left(
|\mathbf g|^2-|\psi|^2
\right)\mathbf I
+
2\operatorname{Re}
\left(
\mathbf g\mathbf g^{\dagger}
\right).
\]

Em uma superfície de raio \(R\), calcule:

\[
\boxed{
\mathbf F
=
-E_0R^2
\int_{4\pi}
\frac{\overline{\mathbf S}}{E_0}
\cdot\widehat{\mathbf e}_r\,d\Omega.
}
\]

O sinal deve ser derivado da orientação da superfície de controle, não
escolhido para forçar concordância.

Escolha:

\[
a<R<\min_{j\ne i}(r_{ij}-a),
\]

para que a superfície encerre apenas a partícula-alvo. Os casos usados no
oráculo devem possuir folga geométrica suficiente.

Verifique:

- pelo menos duas resoluções angulares;
- pelo menos dois raios de controle admissíveis;
- invariância da força com o raio;
- convergência da quadratura;
- concordância com a fórmula analítica para
  \(\mathcal F[a]\), \(\mathcal F[b]\) e \(\mathcal F[c]\).

Use tolerâncias justificadas pelo estudo de refinamento. Como alvo, exija erro
RMS relativo \(\le10^{-5}\) para forças numericamente resolvidas e erro
absoluto normalizado compatível com a precisão da quadratura para forças
nulas. Se esse alvo não for alcançado, refine ou diagnostique; não relaxe
silenciosamente a tolerância.

---

## 6. Testes científicos obrigatórios

Crie, no mínimo:

```text
tests/test_mie_multiparticle.py
tests/test_complete_force.py
tests/test_model_e.py
tests/test_t11_artifacts.py
```

Um auxiliar de oráculo pode ficar em `tests/oracles/` ou em um módulo de
validação claramente não usado pela API de produção.

### 6.1 Sistema global

Teste:

1. \(N=1\): \(b=a\) e \(d=Da\);
2. correspondência material: \(b=a\), \(d=0\) e força de interação nula;
3. identidade direta \(b=a+Ud\);
4. identidade \(d=Db\);
5. resíduo relativo do sistema;
6. equivalência com
   \((I-DU)d=Da\);
7. base planar versus base completa;
8. zeros exatos nos canais inativos;
9. limite rígido \(f_1=1\);
10. permutação de partículas;
11. translação comum no plano \(xy\);
12. rotação comum no plano \(xy\);
13. rejeição de sobreposição e entradas inválidas.

### 6.2 Funcional de força

Teste:

1. saída real e com forma `(3,)`;
2. escala linear em `energy_density`;
3. energia nula;
4. invariância sob fase global:

   \[
   \mathcal F[e^{i\alpha}b]=\mathcal F[b];
   \]

5. coeficientes nulos;
6. índices \(m\) de borda;
7. soma somente até \(M-1\);
8. força externa nula no centro do plano nodal;
9. simetria de um dímero no eixo \(x\);
10. covariância por rotação no plano;
11. concordância com o oráculo de tensor de radiação.

### 6.3 Redução ao observável Rayleigh do Modelo D

Este teste é obrigatório para auditar normalização, fase e conjugação.

Para um caso nodal com o solver \(L=1\) do Modelo D:

1. obtenha o campo incidente externo \(a\);
2. reexpanda o campo espalhado pelas outras partículas até \(n=2\), obtendo
   \(c\);
3. forme \(b=a+c\);
4. use, apenas nesse teste,

   \[
   s_1=i\frac{f_1}{6}(ka)^3,
   \qquad
   s_0=s_2=0;
   \]

5. calcule

   \[
   \mathcal F[b]-\mathcal F[a]-\mathcal F[c].
   \]

Essa parte cruzada deve reproduzir a força do Modelo D em \(L=1\), dentro de
tolerância de ponto flutuante.

Não compare a força completa do Modelo E com o Modelo D esperando igualdade:
o Modelo E contém coeficientes exatos, multipolos superiores e
\(\mathbf F_{\mathrm{ss}}\).

### 6.4 Decomposição

Teste, componente por componente:

\[
\mathbf F_{\mathrm{total}}
-
\mathbf F_{\mathrm{external}}
=
\mathbf F_{\mathrm{ext-sc}}
+
\mathbf F_{\mathrm{ss}}.
\]

Teste também:

- caso \(c=0\);
- caso com \(a=0\) sintético, no qual resta somente
  \(\mathbf F_{\mathrm{ss}}\);
- invariância da decomposição sob fase global;
- ausência de `NaN` e `inf`;
- tratamento explícito de razões não aplicáveis quando a força de referência
  é numericamente nula.

### 6.5 Simetrias físicas

Para duas esferas idênticas em
\((\pm d/2,0,0)\):

\[
F_{1x}^{\mathrm{int}}=-F_{2x}^{\mathrm{int}},
\qquad
F_{iy}^{\mathrm{int}}=F_{iz}^{\mathrm{int}}=0
\]

dentro da tolerância numérica.

Para configurações rotacionadas ou permutadas, compare vetores assinados, não
somente magnitudes.

Não imponha uma lei de ação e reação em clusters genéricos: o campo externo
pode trocar momento com o conjunto.

---

## 7. Convergência e pequena campanha de validação

Crie:

```text
scripts/analyze_t11_model_e.py
```

O script deve ser determinístico e produzir somente uma campanha compacta de
validação. Não reutilize o holdout da T08 como calibração e não inicie a T12.

### 7.1 Casos mínimos

Inclua aproximadamente seis casos:

1. dímero no eixo \(x\), \(d/a=2.5\), \(ka=0.1\), \(f_0=0\), \(f_1=0.8\);
2. dímero diagonal, \(d/a=4\), \(ka=0.05\), \(f_1=0.4\);
3. dímero rígido, \(d/a=3\), \(ka=0.1\), \(f_1=1\);
4. trímero equilátero com lado \(3a\);
5. trímero escaleno não sobreposto;
6. quarteto irregular não sobreposto.

Use \(f_0=0\) nos casos principais para permanecer no domínio científico
atual. Um caso adicional com \(f_0\ne0\) pode ser usado somente para testar a
API geral.

As coordenadas devem ser registradas em documentação ou no CSV, e todas as
geometrias devem ser verificadas quanto a não sobreposição.

### 7.2 Sequência multipolar

Use inicialmente:

\[
L_{\max}=2,3,4,5,6,7.
\]

Se algum caso não estiver convergido, estenda de forma explícita até no
máximo \(L_{\max}=9\) e registre a extensão. Não rotule um caso não confirmado
como divergente.

Para uma força \(\mathbf F_L\), defina a mudança sucessiva:

\[
u_L
=
\frac{
F_{\mathrm{RMS}}(\mathbf F_L-\mathbf F_{L-1})
}{
\max[
F_{\mathrm{RMS}}(\mathbf F_L),
F_{\mathrm{RMS}}(\mathbf F_{L-1})
]
}.
\]

Quando o denominador for numericamente nulo, use uma flag de aplicabilidade e
uma métrica absoluta normalizada por \(a^2E_0\); não produza `NaN` ou `inf`
sem uma flag correspondente.

Considere uma quantidade confirmada somente quando as duas últimas mudanças
sucessivas aplicáveis satisfizerem:

\[
\boxed{u_L\le10^{-5}.}
\]

Avalie separadamente:

- força total;
- força de interação;
- termo `external-scattered`;
- termo `scattered-scattered`.

Um total convergido não prova que uma diferença pequena esteja resolvida.

### 7.3 Artefatos obrigatórios

Gere:

```text
results/data/t11_model_e_convergence.csv
results/data/t11_force_oracle.csv
results/data/t11_force_decomposition.csv
results/figures/t11_model_e_validation.png
```

O CSV de convergência deve incluir, no mínimo:

- `case_id`;
- geometria e \(N\);
- \(ka\), \(f_0\), \(f_1\), \(d_{\min}/a\);
- `lmax`;
- resíduo do sistema;
- número de condição;
- forças RMS normalizadas;
- mudanças sucessivas;
- flags de convergência e aplicabilidade;
- uso ou não da simetria planar.

O CSV do oráculo deve incluir:

- caso, partícula e componente;
- \(L_{\max}\);
- raio de controle;
- resolução angular;
- força analítica;
- força integrada;
- erro absoluto e relativo;
- flag de força numericamente resolvida.

O CSV de decomposição deve incluir:

- forças total, externa e de interação;
- termos `external-scattered` e `scattered-scattered`;
- resíduo de reconstrução;
- razão
  \(F_{\mathrm{RMS}}(F_{\mathrm{ss}})
  /F_{\mathrm{RMS}}(F_{\mathrm{int}})\)
  apenas quando aplicável.

A figura deve ter, preferencialmente, três painéis:

1. convergência multipolar da força de interação;
2. erro do oráculo por tensor de radiação;
3. amplitudes `external-scattered` e `scattered-scattered` nos casos de
   validação.

Use eixos, unidades, legendas e flags de aplicabilidade legíveis. Inspecione
visualmente a figura antes de concluir.

Execute o script duas vezes e confirme determinismo byte a byte dos quatro
artefatos.

---

## 8. Documentação obrigatória

Crie:

```text
TAREFA_T11_MODELO_E_REFERENCIA_COMPLETA.md
```

O documento deve explicar:

- objetivo e limites do Modelo E;
- equação global \((I-UD)b=a\);
- relação \(d=Db\);
- uso da T-matrix exata da T10;
- fórmula completa de força;
- conversão entre \(E_{\mathrm{LAS}}\) e o \(E_0\) do projeto;
- definição de força externa e de interação;
- decomposição em `external-scattered` e `scattered-scattered`;
- distinção entre \(\mathcal F[c]\) e \(2s_ns_{n+1}^{*}\);
- truncamento multipolar da força;
- oráculo independente por tensor de radiação;
- testes de simetria e redução Rayleigh;
- resultados quantitativos da pequena campanha;
- critérios de convergência;
- limitações remanescentes;
- escopo da futura T12.

Atualize apenas o necessário em:

- `README.md`;
- `TASKS.md`;
- `docs/CONVENTIONS.md`;
- `docs/DECISIONS.md`;
- `docs/HANDOFF.md`;
- `src/acoustic_ms/__init__.py`.

Registre explicitamente:

\[
\boxed{
\text{Modelo E}
=
\text{Mie exato}
+
\text{múltiplo espalhamento global}
+
\text{força multipolar completa}.
}
\]

Também registre:

\[
\boxed{
\text{convergência interna do Modelo E}
\ne
\text{validação dos limiares de }\rho_1.
}
\]

A segunda questão pertence à T12-T14.

---

## 9. Restrições de escopo

Nesta tarefa, não:

- altere nenhuma equação ou API dos Modelos A-D;
- conecte silenciosamente Mie ao `multipolar_solver.py` do Modelo D;
- substitua os coeficientes de Rayleigh do Modelo D;
- reutilize `force.py` ou `model_d.py` como fórmula de força do Modelo E;
- copie a fórmula parcial do Modelo D para o novo modelo;
- modifique `rho_foundation.py`;
- refaça a calibração de \(\rho_1\);
- altere os limiares ou o ajuste empírico da T08;
- modifique qualquer CSV ou figura das T01-T10;
- execute a campanha de 20-30 sentinelas da T12;
- crie ou abra o holdout cego da T13-T14;
- introduza absorção, viscosidade, elasticidade sólida, paredes ou esferas
  não idênticas;
- use Mathematica como dependência;
- use notebooks como implementação principal;
- esconda falhas por clipping, regularização ou tolerâncias não documentadas;
- faça commit ou push.

É permitido reutilizar:

- ordenamento de modos;
- harmônicos esféricos;
- coeficientes de translação e Gaunt já validados;
- validação geométrica;
- coeficientes exatos de Mie da T10;
- métricas RMS documentadas.

Se uma mudança fora do escopo parecer indispensável, pare e explique antes de
realizá-la.

---

## 10. Critérios de go/no-go da T11

A T11 só estará tecnicamente aprovada se:

1. a base inicial for o commit `153403c`;
2. os 258 testes-base passarem antes das alterações;
3. o sistema \((I-UD)b=a\) for implementado sem inversão explícita;
4. \(d=Db\) e \(b=a+Ud\) forem verificados;
5. a forma alternativa para \(d\) concordar com a solução principal;
6. a base reduzida concordar com a base completa;
7. a T-matrix da T10 for reutilizada sem duplicação;
8. a força completa seguir as fórmulas e a normalização deste prompt;
9. a redução Rayleigh reproduzir exatamente o observável \(L=1\) do Modelo D;
10. a integração independente do tensor de radiação concordar com a força
    analítica;
11. a identidade da decomposição fechar numericamente;
12. os casos transparentes, rígidos e de partícula única forem aprovados;
13. simetrias, rotações e permutações forem aprovadas;
14. os casos sentinela atingirem o critério multipolar ou forem explicitamente
    marcados como não confirmados;
15. os quatro artefatos forem determinísticos;
16. a figura for inspecionada;
17. todos os testes antigos e novos passarem com warnings como erros;
18. `git diff --check` não apontar problemas;
19. nenhum artefato anterior tiver sido alterado;
20. a documentação não apresentar o Modelo E como incluindo física fora do
    fluido ideal e das esferas sem perdas.

Se o oráculo de tensor de radiação não concordar com a fórmula analítica,
considere isso um **no-go**. Não prossiga para a T12 e não escolha o sinal ou
prefator que produza o resultado desejado. Audite:

- normalização do potencial;
- conversão de energia;
- conjugação;
- orientação da superfície;
- índices \(m\);
- definição de \(h_n^{(1)}\);
- deslocamento `target <- source`.

---

## 11. Verificação final

Ao final, execute:

```bash
.venv/bin/python scripts/analyze_t11_model_e.py
.venv/bin/python scripts/analyze_t11_model_e.py
.venv/bin/python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
git diff --name-only
```

Calcule e registre os SHA-256 dos quatro novos artefatos. Compare os hashes
dos resultados antigos com o inventário inicial.

Não use `git add .`. Não faça commit nem push. Deixe as alterações prontas
para auditoria.

---

## 12. Relatório final ao usuário

Apresente:

1. resumo do Modelo E implementado;
2. commit-base confirmado;
3. equação global realmente resolvida;
4. fórmulas de força e normalização usadas;
5. lista exata de arquivos criados e modificados;
6. quantidade total de testes e resultado da suíte;
7. maior resíduo do sistema;
8. maior diferença entre base completa e reduzida;
9. maior erro do oráculo de tensor de radiação;
10. erro da redução Rayleigh para o Modelo D;
11. maior resíduo da decomposição;
12. estado de convergência de cada caso da campanha;
13. amplitude relativa do termo `scattered-scattered`, sem generalizar além
    dos casos de validação;
14. hashes dos quatro artefatos;
15. confirmação de que os artefatos das T01-T10 permaneceram inalterados;
16. limitações e eventuais casos não confirmados;
17. estado final do Git;
18. veredito objetivo `go` ou `no-go` para iniciar a T12.

Não faça commit nem push. Deixe a T11 pronta para revisão científica e
computacional.
