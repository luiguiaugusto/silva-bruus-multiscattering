# T05.1c — limpeza documental final

## Objetivo e correções

A T05.1c removeu registros antigos contraditórios ou duplicados, corrigiu fragmentos matemáticos malformados, eliminou o desvio RMS incorreto e consolidou a notação de \(\sqrt{2}\), \(f_1\), \(N=3\) e \(L_{\max}=1\). A especificação T05.1 passou a registrar os cinco hashes oficiais.

## Escopo, verificações e resultado

O escopo é exclusivamente documental: os quatro arquivos de tarefa T05.1, `TASKS.md`, `docs/CONVENTIONS.md`, `docs/DECISIONS.md` e `docs/HANDOFF.md`. Código, testes, scripts, dependências e artefatos são protegidos. O aceite exige 92 testes sem warnings, busca vazia pelos defeitos conhecidos, ausência de caracteres ASCII de controle, hashes intactos, `git diff --check` e diff apenas documental. O resultado preservou todos os arquivos científicos. Mensagem de commit: `docs: finalize T05.1 audit cleanup`. T06 não foi iniciada.
