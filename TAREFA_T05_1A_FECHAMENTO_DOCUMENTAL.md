# T05.1a — fechamento documental da auditoria T05.1

## Motivo e escopo

A T05.1a auditou documentalmente a T05.1 no commit `c9967adfee0c4c2266529679631df8f1c89f5131`. As pendências encontradas eram registros incompletos de métricas, geometrias, artefatos, ambiente e limitações. O escopo foi exclusivamente documental: arquivos de tarefa, `TASKS.md` e documentos em `docs/`; código, testes, scripts, resultados, figuras, dependências e forças foram proibidos.

## Verificações e resultado

Foram verificados 92 testes sem warnings, os hashes protegidos e a ausência de alterações científicas. Registrou-se que a T05.1 alterou apenas métricas derivadas e artefatos correspondentes; A, B e C permaneceram inalterados. A auditoria também exige comparar determinismo binário no mesmo ambiente numérico: versões diferentes podem mudar últimos dígitos ou bytes do PNG sem tornar o método aleatório.

Os hashes preservados são T03 `7e02a41ccf3832d233d0e9720f7567ab4eef72ec680df65070f3a687f23fac6a`, T04 `15ee057e2540e7b5f715fa2da4ba13d7f9ed880e0c48ac3cd341f643a5fa37a5`, regressão T05 `e422fff4b12939cc4ea995f03dd04d90f92611f9539549d93a317a6fedaf4ae1`, sweep `dff96cf80380b373b1e9ceab4ef2533df9814553cd8f4c805e8353de6fea50b1` e figura `5327a95c2ccc00151d4389189905feb4b988ea35d8107585f8b9e262ea460d62`.
