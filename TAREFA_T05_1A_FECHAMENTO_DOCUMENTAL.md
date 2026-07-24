# T05.1a — fechamento documental da auditoria T05.1

## Motivo, escopo e aceite

T05.1a documenta a auditoria aprovada no commit `c9967adfee0c4c2266529679631df8f1c89f5131`. Pendências são somente registros técnicos: tarefas, convenções, decisões e handoff. São permitidos somente estes dois arquivos de tarefa, `TASKS.md`, `docs/CONVENTIONS.md`, `docs/DECISIONS.md` e `docs/HANDOFF.md`; código, testes, scripts, artefatos, dependências, forças e métricas são proibidos. Verificar testes, hashes sem regenerar, `git diff --check` e diff apenas documental. Commit: `docs: close T05.1 audit record`; não iniciar T06.

## Números e reprodutibilidade

A auditoria tem 92/92 testes sem warnings; erro de força do oráculo \(\sim2.72	imes10^{-16}\), \(s_{10}\sim6.97	imes10^{-17}\), modo proibido \(\sim1.04	imes10^{-20}\), resíduo \(\sim4.58	imes10^{-16}\), condicionamento \(\sim1.403\), e desvio RMS \(\sqrt2\) \(4.44	imes10^{-16}\). Hashes oficiais: T03 `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`; T04 `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`; regressão T05 `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`; sweep `dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1`; figura `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`. Mesmo ambiente é byte-idêntico; versões diferentes podem alterar bytes, sem aleatoriedade.
