# T11.1 — estabilização numérica do Modelo E

## 1. Objetivo

A T11.1 corrige exclusivamente o caminho numérico do solver multipolar exato
do Modelo E introduzido na T11. As equações físicas, os coeficientes de Mie,
a translação, a força completa e seus quatro canais permanecem inalterados.

A base auditada é o commit
9b2d07719e4f29c2d07f4073939787554b2e966c. O problema observado era o
condicionamento artificialmente extremo da formulação no campo incidente
efetivo em ordens altas, apesar de o problema físico permanecer bem
condicionado após balanceamento modal.

## 2. Sistemas equivalentes

Com \(a\) representando o campo incidente externo, \(b\) o campo incidente
efetivo, \(d=Db\) o campo espalhado, \(D\) a matriz diagonal de espalhamento e
\(U\) a matriz de translação, as duas formas físicas não balanceadas são

\[
A_b b = a,\qquad A_b=I-UD,
\]

e

\[
A_d d = Da,\qquad A_d=I-DU.
\]

Defina a raiz complexa principal, modo a modo,

\[
S=D^{1/2}.
\]

O sistema de produção da T11.1 é

\[
A_q q=Sa,\qquad A_q=I-SUS,\qquad q=Sb.
\]

Somente numpy.linalg.solve é usado. Não se formam inversas, pseudoinversas,
mínimos quadrados ou divisões por \(S\).

## 3. Reconstrução sem divisão

Depois de resolver o sistema balanceado, os campos físicos são reconstruídos
por

\[
d=S q,
\]

\[
b=a+Ud.
\]

Essa reconstrução evita \(S^{-1}q\), preserva exatamente materiais
transparentes e não elimina modos por limiar de magnitude. A única redução de
base permitida continua sendo a simetria planar exata, \(\ell+m\) ímpar.

## 4. Compatibilidade da API

Os atributos públicos legados mantêm seus significados:

- system_matrix: \(A_b\);
- right_hand_side: \(a\);
- condition_number: \(\kappa(A_b)\);
- residual_relative: resíduo relativo avaliado na equação legada.

Foram acrescentados objetos explícitos para \(A_b\), \(A_d\), \(A_q\), seus
lados direitos, diagonais \(D\) e \(S\), coeficientes balanceados, números de
condição e diagnósticos. O campo textual production_solver vale
balanced_sqrt.

## 5. Diagnósticos

O erro retroativo balanceado é

\[
\eta_q=
\frac{\lVert A_q q-Sa\rVert}
{\lVert A_q\rVert\lVert q\rVert+\lVert Sa\rVert}.
\]

Os fechamentos físicos são

\[
r_b=
\frac{\lVert b-a-Ud\rVert}
{\lVert b\rVert+\lVert a\rVert+\lVert Ud\rVert},
\]

\[
r_d=
\frac{\lVert d-Db\rVert}
{\lVert d\rVert+\lVert Db\rVert}.
\]

Quando o denominador é exatamente zero, o diagnóstico retorna zero. Nenhum
piso absoluto é empregado.

## 6. Auditoria das três formulações

A campanha conserva os seis casos da T11 e usa exatamente
\(L_{\max}=2,\ldots,9\), totalizando 48 linhas. A forma legada é resolvida
diretamente apenas para auditoria. A forma espalhada é resolvida diretamente e
recebe uma etapa de refinamento residual. O caminho de produção resolve apenas
\(A_q\).

São comparados \(b\), \(d\) e os quatro canais de força: total, interação,
external–scattered e scattered–scattered. Grandezas resolvidas usam erro
relativo; grandezas numericamente não resolvidas usam diferença absoluta e uma
coluna de aplicabilidade.

## 7. Oráculo de alta precisão

Dois sentinelas em \(L_{\max}=9\), dimer_axis e trimer_scalene, usam mpmath com
80 dígitos decimais. A matriz \(A_q\) e o lado direito oficiais são construídos
em complex128 e convertidos elemento a elemento. O oráculo audita somente a
solução do sistema linear; não constitui uma implementação de Mie ou tradução
em precisão arbitrária.

Comparam-se \(q\), \(d\), \(b\) e os quatro canais de força. O oráculo não é
repetido em testes caros; a campanha determinística registra seus resultados.

## 8. Critérios de aceite

A tarefa exige:

- todos os testes anteriores preservados;
- \(\kappa(A_q)<10\) nos 48 casos;
- erro retroativo balanceado, \(r_b\) e \(r_d\) abaixo de \(10^{-12}\);
- erros de \(q\), \(d\) e \(b\) contra alta precisão abaixo de \(10^{-11}\);
- canais de força resolvidos contra alta precisão abaixo de \(10^{-10}\);
- forças balanceadas contra a forma espalhada abaixo de \(10^{-9}\);
- oráculo de tensão abaixo de \(10^{-10}\);
- artefatos determinísticos;
- artefatos T01–T10 byte a byte inalterados.

## 9. Artefatos

A T11.1 cria results/data/t11_1_solver_stability.csv,
results/data/t11_1_high_precision_oracle.csv e
results/figures/t11_1_model_e_stability.png.

Também regenera os quatro artefatos T11 com o caminho estabilizado e acrescenta
diagnósticos ao CSV de convergência sem renomear as colunas legadas.

## 10. Testes

Os testes cobrem as três matrizes, a raiz principal, reconstrução física,
compatibilidade dos aliases, material transparente, equivalência das
formulações, condicionamento em \(L_{\max}=9\), oráculo mpmath independente em
baixa ordem, finitude, contagens, limites dos artefatos e inspeção do código
para impedir inv, pinv ou lstsq no solver de produção.

## 11. Arquivos protegidos e escopo

Nenhum solver A–D, coeficiente de Mie, tradução, força completa, geometria ou
artefato T01–T10 pode mudar. A T11.1 não executa sentinelas T12, não recalibra
\(\rho_1\) e não abre o holdout T13–T14.

## 12. Verificações finais

Os comandos obrigatórios são:

    .venv/bin/python scripts/analyze_t11_model_e.py
    .venv/bin/python scripts/analyze_t11_1_model_e_stability.py
    .venv/bin/python scripts/analyze_t11_model_e.py
    .venv/bin/python scripts/analyze_t11_1_model_e_stability.py
    .venv/bin/python -m pytest -q -W error
    git diff --check
    git status --short
    git diff --stat
    git diff --name-only

Duas execuções consecutivas no mesmo ambiente devem produzir hashes
idênticos. Esta tarefa termina sem commit e sem push.
