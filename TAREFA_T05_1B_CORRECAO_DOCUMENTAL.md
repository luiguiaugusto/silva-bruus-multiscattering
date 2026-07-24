# T05.1b — correção final da documentação T05.1

## Problemas e correções

Após a T05.1a permaneceram especificações excessivamente curtas e caracteres de controle que corrompiam comandos LaTeX. Esta tarefa substitui os registros por texto autossuficiente, restaura matemática com \(\cdots\) e \[\cdots\], registra o desvio RMS correto \(4.44\times10^{-16}\), e completa o handoff sem alterar ciência.

## Escopo, verificação e aceite

São permitidos somente os dois arquivos de T05.1, este arquivo, `TASKS.md`, `docs/CONVENTIONS.md`, `docs/DECISIONS.md` e `docs/HANDOFF.md`. São protegidos `src/`, `tests/`, `scripts/`, `results/`, `README.md` e dependências. O aceite requer 92 testes sem warnings, hashes oficiais intactos, diff apenas documental, ausência de caracteres ASCII de controle e `git diff --check`. Mensagem de commit: `docs: complete T05.1 audit documentation`.
