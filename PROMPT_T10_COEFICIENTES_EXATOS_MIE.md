# Prompt para execução da T10 — Coeficientes exatos de Mie

Você está trabalhando no repositório:

```text
silva-bruus-multiscattering
```

Execute a **T10 — implementação e validação dos coeficientes exatos de espalhamento acústico de uma esfera compressível**.

Esta tarefa começa após o fechamento da T09. A base esperada da `main` é:

```text
729a2d5 feat: add analytical foundation for rho1
```

O objetivo científico da T10 é substituir, em uma implementação isolada de resposta de partícula única,

\[
s_\ell^{\mathrm{Rayleigh}}
\longrightarrow
s_\ell^{\mathrm{Mie}},
\]

e quantificar o erro causado pelos coeficientes assintóticos de Rayleigh no domínio atual do projeto. **Não integre ainda os coeficientes de Mie ao Modelo D e não implemente o Modelo E.**

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
   - `src/acoustic_ms/scattering.py`;
   - `src/acoustic_ms/multipolar_scattering.py`;
   - `src/acoustic_ms/multipolar_solver.py`;
   - os testes associados a esses módulos.

2. Consulte as referências locais:

   - `papers/Acoustic_Interaction_Force (1)(2).pdf`, especialmente Eq. (10) e Apêndice A;
   - `papers/2014-silva-bruus(2).pdf`, especialmente as definições de \(f_0\) e \(f_1\).

3. Execute:

   ```bash
   git status --short
   git rev-parse --short HEAD
   git log -1 --oneline
   ```

4. O diretório de trabalho deve estar limpo e o `HEAD` deve ser `729a2d5`. Se estiver apenas desatualizado e sem mudanças locais, use:

   ```bash
   git pull --ff-only origin main
   ```

5. Se houver mudanças locais, conflitos ou uma base diferente após o `pull`, **não descarte nem sobrescreva nada**. Pare e informe o estado encontrado.

6. Registre os hashes dos artefatos já versionados sob `results/data/` e `results/figures/`, para confirmar ao final que nenhum resultado anterior foi alterado.

7. Instale o projeto, se necessário, e execute a suíte-base:

   ```bash
   .venv/bin/python -m pip install -e ".[dev,plot]"
   .venv/bin/python -m pytest -q -W error
   ```

O resultado-base esperado é:

```text
217 passed
```

Se a suíte-base falhar, não comece a implementação antes de diagnosticar e relatar a falha.

---

## 2. Escopo científico exato

Implemente os coeficientes parciais exatos para uma esfera fluida, compressível, homogênea e sem perdas, imersa em um fluido ideal ilimitado.

Mantenha:

- a convenção temporal \(e^{-i\omega t}\);
- harmônicos esféricos com a convenção já documentada;
- Hankel esférica de primeira espécie;
- a definição do campo espalhado usada no projeto;
- coeficientes \(s_\ell\) diagonais e independentes de \(m\) para uma esfera;
- propriedades reais e positivas, exceto pelo limite rígido tratado separadamente.

Defina:

\[
x=ka,
\qquad
\widetilde\rho=\frac{\rho_p}{\rho_0},
\qquad
\widetilde\kappa=\frac{\kappa_p}{\kappa_0},
\]

\[
y=k_pa=x\sqrt{\widetilde\rho\,\widetilde\kappa},
\qquad
\beta=
\frac{\rho_0k_p}{\rho_pk}
=
\sqrt{\frac{\widetilde\kappa}{\widetilde\rho}}.
\]

Com

\[
h_\ell^{(1)}(x)=j_\ell(x)+i\,y_\ell(x),
\]

implemente

\[
\boxed{
s_\ell^{\mathrm{Mie}}
=
-
\frac{
\beta j_\ell(x)j_\ell'(y)
-j_\ell(y)j_\ell'(x)
}{
\beta h_\ell^{(1)}(x)j_\ell'(y)
-j_\ell(y){h_\ell^{(1)}}'(x)
}
}
\]

para \(0\le\ell\le L_{\max}\). Todas as linhas indicadas por prima são derivadas em relação ao argumento completo da função.

Use as relações:

\[
\widetilde\kappa=1-f_0,
\]

\[
\widetilde\rho=
\frac{2+f_1}{2(1-f_1)}
\quad\text{para}\quad -2<f_1<1,
\]

\[
\frac{c_p}{c_0}
=
\frac{1}{\sqrt{\widetilde\rho\,\widetilde\kappa}}.
\]

O caso \(f_1=1\) corresponde a \(\widetilde\rho\to\infty\). Não use uma densidade artificialmente grande. Implemente diretamente o limite de esfera rígida:

\[
\boxed{
s_\ell^{\mathrm{rigid}}
=
-
\frac{j_\ell'(x)}
{{h_\ell^{(1)}}'(x)}
}.
\]

Não faça clipping de \(f_1\) nem substitua silenciosamente valores próximos de 1 pelo limite rígido. Apenas \(f_1=1\) deve selecionar explicitamente esse limite.

---

## 3. APIs requeridas

Crie um módulo científico novo, preferencialmente:

```text
src/acoustic_ms/mie_scattering.py
```

Forneça APIs públicas, documentadas e testadas, equivalentes a:

```python
material_ratios_from_contrasts(f0, f1)

fluid_sphere_mie_scattering_coefficients(
    ka,
    density_ratio,
    compressibility_ratio,
    lmax,
)

rigid_sphere_scattering_coefficients(
    ka,
    lmax,
)

mie_scattering_coefficients_from_contrasts(
    ka,
    f0,
    f1,
    lmax,
)
```

Os nomes podem ser refinados se houver uma justificativa clara, mas mantenha separadas:

1. a conversão entre contrastes e propriedades materiais;
2. a solução exata da esfera fluida;
3. o limite rígido;
4. o wrapper em termos de \(f_0,f_1\).

Requisitos:

- retornar um `numpy.ndarray` complexo indexado por \(\ell=0,\ldots,L_{\max}\);
- permitir `lmax=0`;
- validar escalares reais, finitos e domínios físicos;
- exigir \(ka>0\);
- exigir \(\widetilde\rho>0\) e \(\widetilde\kappa>0\) para a esfera fluida;
- no wrapper de contrastes, exigir \(f_0<1\) e \(-2<f_1\le1\);
- rejeitar booleanos onde inteiros ou escalares físicos são esperados;
- retornar exatamente zero no caso de correspondência material exata \(\widetilde\rho=\widetilde\kappa=1\), evitando ruído de cancelamento;
- usar `scipy.special.spherical_jn` e `spherical_yn`, inclusive suas derivadas;
- não adicionar dependências novas sem necessidade demonstrada.

A rotina de partícula única pode aceitar qualquer \(ka>0\) numericamente suportado. Entretanto, toda análise de produção da T10 deve permanecer no domínio do projeto,

\[
10^{-3}\le ka\le0{,}1.
\]

Documente que escolher \(L_{\max}\) suficiente continua sendo responsabilidade do usuário fora desse domínio.

Exporte apenas as APIs públicas necessárias em `src/acoustic_ms/__init__.py`.

---

## 4. Limites de Rayleigh que devem ser recuperados

Compare a nova implementação com a rotina existente:

```python
rayleigh_multipolar_scattering_coefficients
```

Não altere essa rotina.

Os limites obrigatórios são:

\[
s_0^{\mathrm R}
=
-i\frac{f_0}{3}(ka)^3,
\]

e, para \(\ell\ge1\),

\[
s_\ell^{\mathrm R}
=
i
\frac{
3\ell f_1
}{
(2\ell-1)!!(2\ell+1)!!
\left[2(2\ell+1)-(\ell-1)f_1\right]
}
(ka)^{2\ell+1}.
\]

Verifique explicitamente as expressões do Apêndice A da referência de 2026:

\[
s_1=
i\frac{f_1}{6}(ka)^3+O((ka)^5),
\]

\[
s_3=
i\frac{f_1}{350(7-f_1)}(ka)^7+O((ka)^9),
\]

\[
s_5=
i\frac{f_1}{1309770(11-2f_1)}(ka)^{11}
+O((ka)^{13}).
\]

O objetivo não é verificar apenas um ponto: demonstre numericamente a convergência assintótica quando \(ka\to0\), incluindo a ordem esperada da primeira correção relativa.

---

## 5. Validação independente obrigatória

Não valide a fórmula de determinantes apenas copiando a mesma expressão em um teste.

Crie, exclusivamente nos testes, um oráculo independente que resolva para cada \(\ell\) o sistema linear \(2\times2\) obtido das condições de contorno em \(r=a\):

1. continuidade da pressão;
2. continuidade da velocidade radial.

O oráculo deve resolver simultaneamente o coeficiente espalhado \(s_\ell\) e a amplitude interna, usando `numpy.linalg.solve`, e deve ser comparado à implementação de produção em pontos representativos.

Inclua, no mínimo, os seguintes grupos de testes:

### 5.1 Correspondência material

\[
\widetilde\rho=\widetilde\kappa=1
\quad\Rightarrow\quad
s_\ell=0
\]

para todas as ordens testadas.

### 5.2 Condições de contorno

Para diferentes materiais, valores de \(ka\) e ordens \(\ell\), verifique separadamente que os resíduos de:

- pressão;
- velocidade radial;

estão próximos da precisão de máquina.

### 5.3 Oráculo \(2\times2\)

Compare todos os coeficientes de produção ao oráculo independente para pelo menos:

- um material com \(f_0\ne0\) e \(f_1\ne0\);
- um caso com compressibilidade casada, \(f_0=0\);
- \(ka=0{,}01\), \(0{,}05\) e \(0{,}1\);
- \(\ell=0,\ldots,5\).

### 5.4 Limite de Rayleigh

Demonstre que

\[
\frac{s_\ell^{\mathrm{Mie}}}{s_\ell^{\mathrm R}}\to1
\]

para coeficientes Rayleigh não nulos quando \(ka\to0\). Verifique também as potências \(2\ell+1\) e os casos \(\ell=1,3,5\) publicados no Apêndice A.

### 5.5 Limite rígido

Verifique:

- a fórmula direta \(-j_\ell'/h_\ell'\);
- a aproximação por uma sequência de densidades crescentes, apenas como teste de limite;
- o limite Rayleigh \(f_1\to1\) nas ordens positivas.

Evite escolher pontos internos exatamente sobre zeros de funções esféricas na verificação por densidades crescentes.

### 5.6 Unitariedade sem perdas

Para propriedades materiais reais e positivas, teste a identidade compatível com a convenção adotada:

\[
\operatorname{Re}(s_\ell)+|s_\ell|^2=0,
\]

dentro da tolerância numérica apropriada.

### 5.7 Validação de entradas

Teste:

- `NaN` e infinitos;
- \(ka\le0\);
- razões materiais não positivas;
- \(f_0\ge1\);
- \(f_1\le-2\) e \(f_1>1\);
- `lmax` negativo, não inteiro ou booleano;
- entradas complexas onde são exigidos escalares reais.

Use tolerâncias justificadas. Não enfraqueça testes apenas para fazê-los passar.

---

## 6. Campanha numérica da T10

Crie:

```text
scripts/analyze_t10_mie_rayleigh.py
```

A campanha principal deve usar:

\[
f_0=0,
\qquad
f_1\in\{0{,}1,\;0{,}4,\;0{,}8,\;1{,}0\},
\]

\[
ka\in[10^{-3},10^{-1}]
\]

em uma grade logarítmica determinística com pelo menos 101 pontos, e

\[
\ell=0,1,\ldots,5.
\]

Para cada caso, registre, quando matematicamente definido:

- \(ka\);
- \(f_0\) e \(f_1\);
- \(\widetilde\rho\), \(\widetilde\kappa\) e \(c_p/c_0\);
- identificação explícita do limite rígido;
- \(\ell\);
- partes real e imaginária de \(s_\ell^{\mathrm{Mie}}\);
- partes real e imaginária de \(s_\ell^{\mathrm R}\);
- \(|s_\ell^{\mathrm{Mie}}|\) e \(|s_\ell^{\mathrm R}|\);
- erro complexo relativo

  \[
  \varepsilon_{s_\ell}
  =
  \frac{|s_\ell^{\mathrm{Mie}}-s_\ell^{\mathrm R}|}
  {|s_\ell^{\mathrm{Mie}}|};
  \]

- erro relativo de magnitude;
- diferença de fase, somente quando ambos os coeficientes forem não nulos;
- erro absoluto.

Não substitua valores indefinidos por zero. Use `NaN` e uma coluna booleana de aplicabilidade.

Observe que, com \(f_0=0\), o termo Rayleigh \(s_0^{\mathrm R}\) é zero, enquanto o coeficiente exato pode possuir correções de ordem superior. Não use o erro relativo desse canal isoladamente para sugerir relevância física. Registre também o erro absoluto e explique que \(\ell=0\) é inativo pela simetria no problema nodal atual.

Gere:

```text
results/data/t10_mie_rayleigh_validation.csv
results/data/t10_mie_rayleigh_summary.csv
results/figures/t10_mie_rayleigh_error.png
```

O CSV de resumo deve incluir, no mínimo:

- erro máximo no intervalo para cada par \((f_1,\ell)\);
- erro em \(ka=0{,}1\);
- inclinação logarítmica assintótica observada;
- máximo resíduo das condições de contorno;
- máximo defeito de unitariedade;
- separação entre casos fluidos e limite rígido.

A figura deve ser legível e cientificamente enxuta. Use preferencialmente:

1. um painel com o erro do dipolo dominante \(\ell=1\) em função de \(ka\), para os quatro valores de \(f_1\);
2. um painel que resuma, em \(ka=0{,}1\), os erros das ordens \(\ell=1,\ldots,5\).

Não inclua \(\ell=0\) no painel relativo principal quando o denominador Rayleigh for nulo. Use escalas logarítmicas somente onde os dados forem estritamente positivos.

Os artefatos devem ser determinísticos. Execute o script duas vezes e compare os hashes.

---

## 7. Documentação

Crie:

```text
TAREFA_T10_COEFICIENTES_EXATOS_MIE.md
```

Documente:

- problema científico;
- derivação da fórmula implementada;
- relação entre propriedades materiais e \(f_0,f_1\);
- convenções de fase e Hankel;
- tratamento do limite rígido;
- estratégia de estabilidade numérica;
- oráculo independente;
- testes realizados;
- campanha e resultados quantitativos;
- diferença entre exatidão da resposta de uma esfera e exatidão da força coletiva;
- limitações e próximos passos para a T11.

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
\text{coeficientes exatos de Mie}
\ne
\text{força coletiva completa}
}
\]

A T10 valida a T-matrix diagonal de cada esfera. A integração ao sistema global e a implementação da força completa, incluindo termos `scattered–scattered`, pertencem à T11.

Não afirme que a T10 constitui uma validação independente completa de Silva–Bruus ou do Modelo D.

---

## 8. Restrições de escopo

Nesta tarefa, **não**:

- altere `rayleigh_multipolar_scattering_coefficients`;
- altere os resultados ou equações dos Modelos A, B, C ou D;
- conecte silenciosamente Mie ao `multipolar_solver.py`;
- altere `force.py`, `model_d.py`, `rho_foundation.py` ou as fórmulas de força;
- implemente termos `scattered–scattered`;
- implemente o Modelo E;
- refaça ou recalibre \(\rho_1\);
- altere os limiares, ajustes ou arquivos congelados da T08;
- altere os artefatos da T09;
- execute o holdout da T13–T14;
- introduza viscosidade, absorção, elasticidade sólida, paredes ou partículas não esféricas;
- use Mathematica como dependência;
- crie notebooks como implementação principal;
- faça commit ou push.

Se uma mudança fora desse escopo parecer necessária, pare e explique antes de realizá-la.

---

## 9. Critérios de aceite

A T10 só estará concluída se:

1. a implementação reproduzir a Eq. (10) da referência sob as convenções do projeto;
2. o oráculo independente \(2\times2\) concordar com a implementação;
3. as condições de contorno forem satisfeitas numericamente;
4. o caso sem contraste retornar zero;
5. a identidade de unitariedade for verificada;
6. o limite rígido for tratado analiticamente, sem densidade artificial;
7. os limites Rayleigh de \(\ell=0,\ldots,5\) forem recuperados;
8. as fórmulas publicadas para \(\ell=1,3,5\) forem verificadas;
9. os artefatos da campanha forem determinísticos;
10. a figura tiver sido inspecionada visualmente;
11. todos os testes antigos e novos passarem com warnings como erros;
12. `git diff --check` não apontar problemas;
13. nenhum dado ou resultado anterior tiver sido alterado;
14. a documentação distinguir claramente Mie exato, múltiplo espalhamento global e força completa;
15. não houver commit nem push.

Execute ao final:

```bash
.venv/bin/python scripts/analyze_t10_mie_rayleigh.py
.venv/bin/python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
git diff --name-only
```

Compare também os hashes dos artefatos anteriores com o registro inicial.

---

## 10. Relatório final ao usuário

Ao terminar, apresente:

1. resumo do que foi implementado;
2. fórmula e convenções efetivamente usadas;
3. lista exata de arquivos criados e modificados;
4. quantidade total de testes e resultado da suíte;
5. erros máximos do oráculo, condições de contorno e unitariedade;
6. erros Mie–Rayleigh mais importantes, sobretudo para \(\ell=1\) em \(ka=0{,}1\);
7. resultado do limite rígido;
8. hashes dos três novos artefatos;
9. confirmação de que os artefatos anteriores permaneceram inalterados;
10. limitações remanescentes;
11. estado final do Git;
12. avaliação objetiva sobre a prontidão da T10 para auditoria e commit.

Não faça commit nem push. Deixe a T10 pronta para revisão.
