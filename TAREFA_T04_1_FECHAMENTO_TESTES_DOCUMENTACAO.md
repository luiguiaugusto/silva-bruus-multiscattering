# T04.1 — Fechamento da cobertura de testes e da documentação da T04

## Instrução de uso

Copie esta especificação integralmente para o Codex CLI aberto na raiz do
repositório `silva-bruus-multiscattering`.

A implementação científica da T04 já foi auditada e aprovada. A T04.1 é uma
tarefa de fechamento formal: deve acrescentar a cobertura de testes que faltou
e corrigir a documentação, sem modificar o cálculo da força, o solver
multipolar ou qualquer resultado científico aprovado.

---

## 1. Estado inicial e preparação obrigatória

A T04 aprovada está na `main`, após os commits:

```text
93686e43772cbf46de2f8b42b43479d639a732b6
ad797f58e79e0b3fb91a02d729cf167c4385f47e
```

Antes de alterar qualquer arquivo:

1. execute:

   ```bash
   git status --short
   git branch --show-current
   git rev-parse HEAD
   git log -5 --oneline
   ```

2. preserve qualquer alteração preexistente do usuário;
3. confirme que a branch é a mesma `main` usada na T04 e que o histórico contém
   os dois commits acima;
4. leia integralmente:
   - `AGENTS.md`;
   - `TAREFA_T04_FORCA_INTERACAO_MS_RAYLEIGH.md`;
   - `src/acoustic_ms/force.py`;
   - `src/acoustic_ms/solver.py`;
   - `src/acoustic_ms/__init__.py`;
   - `tests/test_force.py`;
   - `docs/CONVENTIONS.md`;
   - `docs/DECISIONS.md`;
   - `docs/HANDOFF.md`;
   - `TASKS.md`;
5. instale o projeto no ambiente já usado pelo repositório e execute:

   ```bash
   python -m pytest -q -W error
   ```

6. confirme a linha de base de:

   ```text
   75 passed
   ```

Se a árvore estiver suja, a branch estiver incorreta, os commits da T04 não
estiverem no histórico ou a linha de base não tiver 75 testes aprovados, pare e
relate o problema. Não use `git reset --hard`, `git checkout --`, `git clean`,
`--force` ou outro comando destrutivo.

---

## 2. Escopo estrito

A T04.1 deve alterar somente:

- testes da força;
- documentação do estado já implementado;
- o docstring de abertura do pacote;
- o registro da tarefa em `TASKS.md`.

Não altere:

- `src/acoustic_ms/force.py`;
- `src/acoustic_ms/solver.py`;
- os módulos de espalhamento, translação, funções especiais ou Gaunt;
- as fórmulas das T01 e T02;
- os CSVs validados;
- tolerâncias de testes anteriores.

Se algum novo teste revelar um defeito científico real e reproduzível, não
corrija silenciosamente o núcleo. Pare, apresente o menor contraexemplo e
explique o impacto antes de modificar qualquer módulo científico.

---

## 3. Testes faltantes

Faça as adições em `tests/test_force.py`. É permitido reorganizar somente o
código de teste para evitar repetição, desde que os testes já existentes
continuem equivalentes.

### A. Independência da ordem das partículas

Adicione um teste explícito para um par oblíquo, sem introduzir resultados para
três ou mais partículas.

Procedimento:

1. construa um par não alinhado aos eixos;
2. calcule o resultado original;
3. permute a entrada com uma ordem como:

   ```python
   order = np.array([1, 0])
   ```

4. recalcule o resultado;
5. verifique, após aplicar a mesma permutação aos resultados de referência:
   - `forces_xy`;
   - `local_scattered_coefficients`;
   - `solution.coefficients`.

Use `np.allclose` com tolerâncias compatíveis com os testes científicos já
aprovados, por exemplo `rtol=2e-12` e `atol=2e-14`. Não compare apenas a soma
das forças: o teste deve demonstrar a correspondência partícula a partícula.

### B. Modos locais proibidos pela simetria

No par alinhado ao eixo \(x\), o campo espalhado nodal reexpandido até
\(\ell=2\) pode conter, nesta configuração:

\[
(1,0),\qquad(2,-1),\qquad(2,1).
\]

Os modos proibidos são:

\[
(0,0),\quad(1,-1),\quad(1,1),\quad
(2,-2),\quad(2,0),\quad(2,2).
\]

Adicione uma verificação explícita de que esses seis coeficientes são nulos
dentro da precisão numérica para as duas partículas.

Construa os índices com `mode_index(ell, m)`; não use números mágicos. Uma
tolerância absoluta da ordem de \(2\times10^{-13}\) é adequada. Preserve também
o teste já existente:

\[
b_{2,-1}=-b_{2,1}.
\]

Esse é um teste da reexpansão local usada no observável de força. Ele não
promove o solver de espalhamento para \(L_{\max}=2\).

### C. Energia acústica não finita

Adicione testes explícitos para:

```python
energy_density = np.nan
energy_density = np.inf
energy_density = -np.inf
```

Todos devem produzir `ValueError`. Mantenha o teste já existente para energia
negativa. Não suprima warnings nem permita que a falha ocorra posteriormente
por propagação de `NaN`.

### D. Limite inferior inválido de \(f_1\)

A cobertura atual testa \(f_1>1\), mas não testa o outro lado do domínio.
Adicione um caso explícito, por exemplo:

```python
f1 = -2.01
```

e confirme `ValueError`.

Não altere o domínio aprovado:

\[
-2\leq f_1\leq1.
\]

### E. `lmax != 1`

Adicione teste direto da API pública
`solve_rayleigh_nodal_interaction_forces` para valores não suportados, por
exemplo:

```python
lmax = 0
lmax = 2
```

Ambos devem produzir `ValueError`. Isso deve confirmar que:

\[
L_{\max}^{\mathrm{scatter}}=1
\]

continua sendo uma restrição do solver, embora a avaliação local da força use
\(\ell\leq2\).

### F. Formato inválido de `positions_xyz`

Adicione ao menos um teste direto da API de força com formato incompatível com
\((N,3)\), por exemplo:

```python
positions_xyz = [0.0, 0.0, 0.0]
```

ou:

```python
positions_xyz = [[0.0, 0.0]]
```

e confirme `ValueError`. O teste deve passar pela API pública da T04, não chamar
apenas a função privada de validação do solver.

### G. Regras para os novos testes

- Não use `skip`, `xfail` ou supressão de warnings.
- Não altere tolerâncias antigas para esconder falhas.
- Não construa valores esperados chamando a própria função sob teste.
- Não adicione dependências.
- Não crie resultados científicos para \(N\geq3\).
- Preserve todos os 75 testes da linha de base.

A contagem final será maior que 75 e dependerá de como os casos parametrizados
forem coletados. Registre no `HANDOFF.md` a contagem efetivamente informada pelo
`pytest`; não escreva antecipadamente um total estimado.

---

## 4. Correções de documentação

### A. `docs/CONVENTIONS.md`

Na seção da T04, corrija as duas fórmulas cartesianas para preservar
explicitamente a conjugação complexa:

\[
F_x=
\frac{\sqrt{30\pi}}{15}ka^3E_0
\operatorname{Re}
\left[
f_1^*
\left(b_{2,-1}-b_{2,1}\right)
\right],
\]

\[
F_y=
\frac{\sqrt{30\pi}}{15}ka^3E_0
\operatorname{Re}
\left[
-i f_1^*
\left(b_{2,1}+b_{2,-1}\right)
\right].
\]

O código atual já usa `np.conj(f1)` e está correto. Esta alteração é apenas
documental. Registre, se útil, que a API atual restringe \(f_1\) a um escalar
real, mas mantém a conjugação da derivação.

Não mude sinal, prefator, direção, normalização de energia ou interpretação da
força.

### B. `src/acoustic_ms/__init__.py`

Atualize somente o docstring inicial para incluir:

- T01: força de pares original de Silva–Bruus;
- T02: benchmark analítico corrigido de duas partículas;
- T03: solver acoplado Rayleigh com \(L_{\max}=1\);
- T04: força de interação nodal do Modelo C, com
  \(L_{\max}^{\mathrm{scatter}}=1\) e avaliação local até \(\ell=2\).

Não altere `__all__`, nomes públicos, imports ou comportamento da API.

### C. `docs/HANDOFF.md`

Corrija e complete a seção T04.

1. Substitua a afirmação incorreta de `70 tests` pela contagem final realmente
   obtida após a T04.1.
2. Crie uma subseção `T04.1 coverage closure` ou equivalente.
3. Registre:
   - os novos testes de permutação;
   - os seis modos locais proibidos;
   - energia não finita;
   - \(f_1<-2\);
   - `lmax != 1`;
   - formato inválido de `positions_xyz`;
   - as tolerâncias efetivamente usadas.
4. Registre os comandos finais executados.
5. Registre os hashes SHA-256 dos CSVs da T03 e da T04.
6. Confirme que ambos são regenerados deterministicamente e que a T04.1 não
   altera seus bytes.
7. Registre os erros efetivamente recalculados para:
   - oráculo cartesiano independente;
   - referência escalar do par;
   - limite \(f_1\rightarrow0\);
   - limite de grande separação.
8. Registre também:
   - resíduo relativo máximo;
   - maior número de condição;
   - violação máxima de ação–reação.
9. Preserve os quatro valores de regressão no contato.
10. Mantenha explícito que:
    - a T04 retorna força de interação nodal, não força total irrestrita;
    - o solver permanece em \(L_{\max}=1\);
    - \(\ell=2\) é somente ordem de avaliação local;
    - nenhum resultado para \(N\geq3\) foi produzido.

Como verificações de sanidade, a auditoria da T04 obteve valores da ordem de:

| Métrica | Valor auditado |
|---|---:|
| Erro máximo contra a referência escalar | \(4.276\times10^{-16}\) |
| Erro máximo contra o oráculo cartesiano | \(4.924\times10^{-16}\) |
| Erro no limite \(f_1\rightarrow0\) | \(6.374\times10^{-9}\) |
| Erro no limite de grande separação | \(3.891\times10^{-8}\) |
| Resíduo relativo máximo | \(1.498\times10^{-16}\) |
| Maior número de condição | \(1.286\) |
| Violação máxima de ação–reação | \(1.481\times10^{-16}\) |

Recalcule e registre os valores reais produzidos no ambiente atual. Pequenas
diferenças nos últimos algarismos são aceitáveis; não copie os números da
tabela sem verificá-los.

O hash auditado do CSV da T04 é:

```text
15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5
```

Confirme-o após duas regenerações. Para o CSV da T03, registre o hash encontrado
antes e confirme que o mesmo valor permanece depois da regeneração.

### D. `TASKS.md`

Adicione e marque como concluída uma entrada curta para a T04.1, descrevendo-a
como fechamento da cobertura de testes e da documentação. Não antecipe T05 nem
marque tarefas científicas posteriores como concluídas.

---

## 5. Verificação de não regressão científica

Confirme que nenhum módulo científico foi alterado em relação ao commit-base da
T04. Depois das edições, um comando equivalente a:

```bash
git diff --exit-code ad797f58e79e0b3fb91a02d729cf167c4385f47e \
  -- src/acoustic_ms/contrasts.py \
     src/acoustic_ms/silva_bruus.py \
     src/acoustic_ms/corrected_pair.py \
     src/acoustic_ms/multipoles.py \
     src/acoustic_ms/special.py \
     src/acoustic_ms/gaunt.py \
     src/acoustic_ms/translation.py \
     src/acoustic_ms/scattering.py \
     src/acoustic_ms/incident.py \
     src/acoustic_ms/solver.py \
     src/acoustic_ms/force.py
```

deve terminar com código zero e sem diff.

O docstring de `src/acoustic_ms/__init__.py` pode mudar, mas sua API, imports e
`__all__` devem permanecer idênticos.

---

## 6. Verificação final

Execute, nesta ordem:

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q

sha256sum results/data/t03_solver_validation.csv
sha256sum results/data/t04_pair_force_validation.csv

python scripts/validate_t03_solver.py
python scripts/validate_t04_force.py
sha256sum results/data/t03_solver_validation.csv
sha256sum results/data/t04_pair_force_validation.csv

python scripts/validate_t04_force.py
sha256sum results/data/t04_pair_force_validation.csv

python -m pytest -q -W error
git diff --check
git status --short
git diff --stat
```

Confirme que:

- os 75 testes anteriores continuam passando;
- todos os novos testes passam;
- não existem warnings;
- os hashes são estáveis;
- o CSV da T04 conserva o hash auditado;
- o CSV da T03 não muda;
- nenhum módulo científico foi modificado;
- `src/acoustic_ms/__init__.py` mudou apenas no docstring;
- não foi produzido resultado para \(N\geq3\).

---

## 7. Entrega, commit e push

Ao finalizar, apresente:

1. o commit-base encontrado;
2. a saída completa de `python -m pytest -q -W error`;
3. a contagem de testes da linha de base, novos e total;
4. os erros numéricos recalculados;
5. os hashes dos CSVs da T03 e T04;
6. a confirmação de determinismo byte a byte;
7. os arquivos alterados;
8. `git diff --stat`;
9. a confirmação explícita de que `force.py`, `solver.py` e os demais módulos
   científicos não mudaram;
10. a confirmação de que o solver continua em \(L_{\max}=1\), que \(\ell=2\)
    permanece somente avaliação local e que não existem resultados para
    \(N\geq3\).

Faça um commit com mensagem semelhante a:

```text
test: close T04 validation gaps
```

Depois, faça o push para a mesma branch remota atualmente utilizada.

Se o push exigir trocar o remote, sobrescrever histórico, usar `--force` ou
resolver divergência não relacionada à T04.1, pare e informe o problema.
