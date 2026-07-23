# Revisão 01 do projeto Silva–Bruus–multiscattering

**Data:** 23 de julho de 2026  
**Objetivo:** registrar a auditoria científica inicial, corrigir o plano de execução e fornecer a primeira tarefa pronta para o Codex CLI.

## 1. Veredito

A ideia científica é sólida e tem potencial de publicação, mas a novidade não deve ser apresentada como “uma formulação para \(N\) esferas”. Uma formulação de múltiplo espalhamento e força para um conjunto de \(N\) esferas em fluido ideal já foi apresentada por Lopes, Azarpeyvand e Silva (2016).

A contribuição defensável é:

> quantificar, decompor e prever a quebra da aproximação pairwise de Silva–Bruus em clusters finitos de partículas Rayleigh.

O trabalho deve separar:

1. correções completas de dois corpos;
2. contribuições genuinamente não pairwise;
3. erro devido ao truncamento multipolar;
4. erro devido à aproximação Rayleigh da resposta de uma esfera;
5. erro decorrente da própria fórmula usada para calcular a força a partir do campo acoplado.

Referências centrais:

- Silva e Bruus, *Physical Review E* 90, 063007 (2014): <https://doi.org/10.1103/PhysRevE.90.063007>
- Lopes, Azarpeyvand e Silva, *IEEE TUFFC* 63, 186–196 (2016): <https://doi.org/10.1109/TUFFC.2015.2494693>
- Sepehrirahnama e Lim, *Physical Review E* 102, 043307 (2020): <https://doi.org/10.1103/PhysRevE.102.043307>
- Silva, *Brazilian Journal of Physics* (2026): <https://doi.org/10.1007/s13538-026-02102-x>

## 2. Correções obrigatórias no plano

### 2.1 Erro de sinal na fórmula de Silva–Bruus

O plano anterior registrou o termo \(kd\sin(kd)\) com sinal positivo dentro do primeiro parêntese. A forma correta, usando

\[
E_0=\frac{\rho_0|v_0|^2}{4}
\]

para amplitudes complexas, é

\[
\mathbf F_{i\leftarrow j}^{\mathrm{SB}}
=
4\pi a^2E_0f_1^2
\left(\frac{a}{d}\right)^4
B(kd)\,\widehat{\mathbf d}_{ij},
\]

\[
B(x)=
-\frac32\left[\cos x+x\sin x\right]
+\frac12x^2\cos x,
\]

onde

\[
\mathbf d_{ij}=\mathbf r_i-\mathbf r_j,\qquad
d=\|\mathbf d_{ij}\|.
\]

Portanto, \(\widehat{\mathbf d}_{ij}\) aponta da partícula-fonte \(j\) para a partícula-prova \(i\). Para \(kd\ll1\), \(B(kd)<0\), e a força é atrativa.

Os limites corretos são

\[
\mathbf F_{\mathrm{SB}}
\sim
-6\pi a^2E_0f_1^2
\left(\frac ad\right)^4
\widehat{\mathbf d},
\qquad kd\ll1,
\]

e

\[
\mathbf F_{\mathrm{SB}}
\sim
2\pi a^2E_0f_1^2(ka)^2
\left(\frac ad\right)^2
\cos(kd)\widehat{\mathbf d},
\qquad kd\gg1,\quad a/d\ll1.
\]

### 2.2 Normalização deve ser única

O artigo de 2014 e o trabalho de 2026 empregam convenções de amplitude/energia que podem produzir um fator dois aparente. O projeto adotará:

\[
p(\mathbf r,t)=\operatorname{Re}\{p(\mathbf r)e^{-i\omega t}\},
\qquad
E_0=\frac{\rho_0|v_0|^2}{4}.
\]

Essa convenção deve aparecer no código, nos testes, nos gráficos e em `docs/CONVENTIONS.md`.

### 2.3 \(ka\), \(kd\) e \(a/d\) não são independentes

Deve-se parametrizar cada configuração por \(ka\) e pelas posições adimensionais \(\mathbf r_i/a\). Então,

\[
kd=(ka)\frac da.
\]

Varreduras não podem escolher simultaneamente \(ka\), \(kd\) e \(d/a\) como variáveis independentes.

### 2.4 O vetor desconhecido do sistema global deve ser definido

Adotaremos coeficientes espalhados \(\mathbf s\) como incógnitas:

\[
\mathbf s_i
=
\mathbf T_i
\left(
\mathbf a_i^{\mathrm{ext}}
+\sum_{j\ne i}\mathbf U_{ij}\mathbf s_j
\right),
\]

\[
(\mathbf I-\mathbf T\mathbf U)\mathbf s
=
\mathbf T\mathbf a^{\mathrm{ext}}.
\]

\(\mathbf U_{ij}\) transforma coeficientes outgoing no centro \(j\) em coeficientes regular no centro \(i\). O vetor

\[
\mathbf r_{ij}=\mathbf r_i-\mathbf r_j
\]

e a convenção temporal \(e^{-i\omega t}\) não podem ser trocados durante o projeto.

### 2.5 A “força completa” ainda precisa ser definida

Resolver o campo acoplado não basta. Também é necessário calcular a força com uma expressão compatível com o nível de aproximação desejado.

A equação de Gor’kov aplicada ao campo incidente efetivo é quadrática. Ao escrever

\[
p_i^{\mathrm{in}}=p_i^{\mathrm{ext}}+\sum_{j\ne i}p_j^{\mathrm{sc}},
\]

aparecem:

1. termos externo–externo, que produzem a força primária;
2. termos cruzados externo–espalhado;
3. termos espalhado–espalhado.

A expressão de dois corpos publicada em 2026 retém os termos cruzados usados em sua ordem assintótica. Ela é um benchmark obrigatório, mas não deve ser assumida automaticamente como a referência mais completa para clusters densos.

O projeto deverá implementar e comparar, em fases posteriores:

- força cruzada, para reproduzir o artigo de 2026;
- força de Gor’kov do campo incidente efetivo menos a força primária, incluindo termos espalhado–espalhado;
- fórmula multipolar de força por fluxo de momento/far field de Lopes et al. (2016), que será a referência geral.

### 2.6 \(L_{\max}=5\) é benchmark, não garantia de convergência

O artigo de 2026 fornece um benchmark analítico truncado em \(L_{\max}=5\). Para o novo paper, a referência numérica deverá satisfazer um critério de convergência, por exemplo

\[
\frac{\|\mathbf F^{(L)}-\mathbf F^{(L-2)}\|}
{\max(\|\mathbf F^{(L)}\|,F_{\mathrm{floor}})}
<\tau_L.
\]

Devem ser testados \(L=1,3,5,7,\ldots\), sobretudo próximo ao contato e para alto contraste. Não se deve chamar \(L=5\) de “solução completa” antes dessa análise.

### 2.7 Separar T-matrix exata e expansão Rayleigh

Serão necessários dois modelos de espalhamento de uma esfera:

- coeficientes assintóticos usados no benchmark analítico de 2026;
- coeficientes exatos de uma esfera compressível em fluido ideal.

Isso separa erro de rescattering, erro multipolar e erro da aproximação de pequena partícula.

### 2.8 Explorar a simetria do plano nodal

Para um cluster inteiro contido no plano \(z=0\), o problema é invariante por \(z\mapsto-z\). Na base esférica usual, sobrevivem os canais compatíveis com a paridade

\[
\ell+m\ \text{ímpar}.
\]

Isso reduz \(L_{\max}=5\) de 36 para 18 modos por partícula no problema planar geral. No truncamento \(L_{\max}=1\), resta apenas o canal dipolar vertical \((\ell,m)=(1,0)\), transformando o primeiro modelo multibody em um sistema escalar de dimensão \(N\).

Essa redução é uma oportunidade importante: permite estudar analiticamente a conexão entre Silva–Bruus, a primeira aproximação de Born e a série de Neumann do operador de rescattering antes de implementar toda a maquinaria multipolar.

### 2.9 O raio espectral não deve ser o único critério

\(\rho(\mathbf T\mathbf U)\) é fisicamente relevante, mas pode ser insuficiente para operadores não normais e não incorpora como o campo externo excita os modos.

Também devem ser avaliados:

\[
\|\mathbf T\mathbf U\|,
\qquad
\frac{\|\mathbf T\mathbf U\,\mathbf s^{(0)}\|}
{\|\mathbf s^{(0)}\|},
\qquad
\|(\mathbf I-\mathbf T\mathbf U)^{-1}\|,
\]

além da métrica geométrica proposta,

\[
\Lambda_i=|f_1|\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3.
\]

O objetivo será descobrir qual indicador melhor prediz o erro da força, e não assumir antecipadamente que \(\Lambda_i\) ou \(\rho\) é universal.

### 2.10 Clusters aleatórios exigem estatística definida

Para cada ensemble deverão ser registrados:

- domínio finito;
- fração de área;
- distância mínima;
- algoritmo de geração;
- semente aleatória;
- número de realizações;
- intervalo de confiança;
- classificação de partículas internas e de borda.

### 2.11 \(N=1000\) não é requisito científico

O paper pode ser forte com \(N\le100\). A escala máxima será determinada pela pergunta física e pela convergência estatística. \(N=500\) ou \(1000\) ficará como demonstração opcional após a validação; não deve atrasar a contribuição principal em \(N=3\) e pequenos clusters.

## 3. Roteiro revisado

1. **T01 — Silva–Bruus pairwise e convenções.**
2. **T02 — Fórmula corrigida de dois corpos de 2026 e reprodução da Figura 2.**
3. **T03 — Modelo nodal \(L_{\max}=1\), sistema escalar para \(N\), série de Neumann e contribuição irredutível de três corpos.**
4. **T04 — Funções especiais, índices e coeficientes de esfera.**
5. **T05 — Operador de translação, com testes diretos do teorema de adição e entradas do Apêndice D do artigo de 2026.**
6. **T06 — Solver multipolar geral; comparação do solver completo de \(N=2\) com o sistema reduzido e a fórmula publicada.**
7. **T07 — Fórmula geral da força e comparação entre força cruzada, Gor’kov efetiva e far field.**
8. **T08 — \(N=3\): cadeia, triângulo equilátero e triângulo escaleno.**
9. **T09 — Convergência em \(L_{\max}\), mapas de erro e critérios de validade.**
10. **T10 — Ensembles finitos e escalabilidade somente quando cientificamente necessários.**

## 4. Tarefa pronta para enviar ao Codex CLI

Copie a especificação abaixo integralmente.

---

# T01 — Infraestrutura e força nodal pairwise de Silva–Bruus

## Objetivo

Crie a infraestrutura inicial do repositório Python e implemente somente:

1. fatores de contraste \(f_0\) e \(f_1\);
2. força de interação pairwise de Silva–Bruus para duas esferas idênticas presas em um plano nodal de pressão;
3. testes físicos, algébricos e numéricos;
4. documentação explícita das convenções.

Não implemente múltiplo espalhamento nesta tarefa.

## Convenções obrigatórias

Use:

\[
p(\mathbf r,t)=\operatorname{Re}\{p(\mathbf r)e^{-i\omega t}\},
\qquad
E_0=\frac{\rho_0|v_0|^2}{4}.
\]

As posições passadas à função vetorial pertencem ao plano nodal e serão vetores 2D.

Para a força sobre a partícula \(i\) produzida pela partícula \(j\), defina

\[
\mathbf d_{ij}=\mathbf r_i-\mathbf r_j,\qquad
d=\|\mathbf d_{ij}\|,\qquad
\widehat{\mathbf d}_{ij}=\frac{\mathbf d_{ij}}d.
\]

A equação a implementar é

\[
\mathbf F_{i\leftarrow j}^{\mathrm{SB}}
=
4\pi a^2E_0f_1^2
\left(\frac{a}{d}\right)^4
\left\{
-\frac32[\cos(kd)+kd\sin(kd)]
+\frac12(kd)^2\cos(kd)
\right\}
\widehat{\mathbf d}_{ij}.
\]

O sinal de \(kd\sin(kd)\) dentro do primeiro colchete é **negativo** após a multiplicação por \(-3/2\). Não use a versão antiga com \(+\frac32kd\sin(kd)\).

Os contrastes são

\[
f_0=1-\frac{\kappa_p}{\kappa_0},
\qquad
f_1=\frac{2(\rho_p/\rho_0-1)}{2\rho_p/\rho_0+1}.
\]

## Estrutura mínima

Crie:

```text
pyproject.toml
README.md
AGENTS.md
TASKS.md
docs/CONVENTIONS.md
docs/DECISIONS.md
docs/HANDOFF.md
src/acoustic_ms/__init__.py
src/acoustic_ms/contrasts.py
src/acoustic_ms/silva_bruus.py
tests/test_contrasts.py
tests/test_silva_bruus.py
```

Use layout `src`, Python 3.11 ou superior, NumPy como dependência de execução e pytest como dependência de desenvolvimento. Não adicione SciPy, SymPy, Jupyter, pandas ou Numba nesta tarefa.

## API mínima

Implemente funções equivalentes a:

```python
monopole_contrast(compressibility_ratio)
dipole_contrast(density_ratio)
nodal_pair_force_magnitude(k, radius, distance, energy_density, f1)
nodal_pair_force_on_probe(probe_xy, source_xy, k, radius, energy_density, f1)
nodal_pair_forces(position_1_xy, position_2_xy, k, radius, energy_density, f1)
```

Os nomes podem ser ajustados apenas se houver justificativa clara no `docs/HANDOFF.md`.

Requisitos:

- aceitar escalares NumPy quando natural;
- documentar unidades SI;
- rejeitar valores não finitos;
- exigir \(k>0\), \(a>0\), \(E_0\ge0\) e \(d\ge2a\);
- exigir razão de densidade positiva;
- aceitar \(f_1\) no intervalo físico limite \([-2,1]\);
- a função `nodal_pair_forces` deve retornar forças iguais e opostas;
- não ocultar a normalização dentro de constantes globais.

## Teste independente pelo potencial

Além de testar a fórmula fechada, verifique numericamente que ela coincide com

\[
F(d)=-\frac{dU}{dd},
\]

onde, na convenção de energia adotada,

\[
U(d)=2\pi E_0 k^3a^6f_1^2\frac{n_1(kd)}{kd},
\]

\[
n_1(x)=-\frac{\sin x}{x}-\frac{\cos x}{x^2}.
\]

O teste deve avaliar vários valores de \(kd\), sem reutilizar internamente a função da força para construir o valor esperado.

## Testes obrigatórios

1. \(f_0=0\) para \(\kappa_p/\kappa_0=1\).
2. \(f_1=0\) para \(\rho_p/\rho_0=1\).
3. \(f_1\to1\) para razão de densidade muito grande.
4. \(f_1\to-2\) para razão de densidade positiva tendendo a zero.
5. Identidade entre força fechada e \(-dU/dd\) em valores representativos de \(kd\).
6. Limite
   \[
   F/F_{\mathrm{near}}\to1,\qquad
   F_{\mathrm{near}}=-6\pi a^2E_0f_1^2(a/d)^4,
   \]
   para \(kd\ll1\).
7. Limite
   \[
   F/F_{\mathrm{far}}\to1,\qquad
   F_{\mathrm{far}}=2\pi a^2E_0f_1^2(ka)^2(a/d)^2\cos(kd),
   \]
   para \(kd\gg1\), escolhendo um ponto afastado de zeros de \(\cos(kd)\).
8. Ação–reação.
9. Invariância por translação dentro do plano nodal.
10. Covariância por rotação dentro do plano nodal.
11. Força nula para \(f_1=0\).
12. Força nula para \(E_0=0\).
13. Força atrativa para \(kd\ll1\).
14. Rejeição de sobreposição \(d<2a\).
15. Rejeição de entradas inválidas e posições coincidentes.

Não teste invariância por uma rotação 3D arbitrária nem por translação na direção \(z\): o campo externo seleciona o eixo \(z\) e o plano nodal.

## Documentação

`docs/CONVENTIONS.md` deve registrar:

- convenção temporal;
- definição de \(E_0\);
- direção de \(\widehat{\mathbf d}_{ij}\);
- interpretação de sinal;
- domínio \(d\ge2a\);
- diferença entre pressão nodal e antinodal;
- relação \(kd=(ka)(d/a)\);
- referência às equações de Silva–Bruus e ao artigo de 2026.

`docs/DECISIONS.md` deve registrar que Python é a implementação oficial e que os notebooks não conterão a única versão de nenhuma rotina científica.

## Fora do escopo

Não implemente:

- Bessel, Hankel ou harmônicos esféricos;
- coeficientes de Gaunt;
- matriz de translação;
- T-matrix;
- solver global;
- fórmula corrigida de dois corpos de 2026;
- gráficos;
- notebooks;
- otimizações.

## Verificação

Execute:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

Todos os testes devem passar sem warnings produzidos pelo código do projeto.

## Entrega

Atualize `docs/HANDOFF.md` com:

1. arquivos criados;
2. equações implementadas;
3. convenções adotadas;
4. comandos executados;
5. resultado dos testes;
6. limitações ou dúvidas;
7. resumo de `git diff --stat`.

Pare após concluir a T01. Não avance para a fórmula corrigida de dois corpos nem para múltiplo espalhamento.

---

## 5. Critério para autorizar a T02

A T02 só será iniciada depois que o relatório do Codex demonstrar:

- uso da fórmula com o sinal correto;
- normalização explícita de \(E_0\);
- teste independente por derivada do potencial;
- testes restritos às simetrias realmente válidas no plano nodal;
- suíte pytest integralmente aprovada;
- nenhuma implementação prematura de múltiplo espalhamento.
