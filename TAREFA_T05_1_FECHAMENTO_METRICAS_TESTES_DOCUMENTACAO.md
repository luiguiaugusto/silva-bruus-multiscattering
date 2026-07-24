# T05.1 -- fechamento de métricas, testes, artefatos e documentação

Registro da especificação de encerramento da T05.1.

Escopo: corrigir somente as métricas de nulidade numérica e RMS vetorial, ampliar a cobertura dos trímeros, regenerar os artefatos T05 e atualizar a documentação. A/B/C, equações, solver, forças, oráculo escalar e regressões científicas permanecem protegidos.

Métricas: `F_scale=max_i(|F_ref,i|,|F_mod,i|)`, `F_tol=128*eps*F_scale`, sem piso absoluto; ângulos de vetores numericamente nulos são NaN; `F_RMS=sqrt(mean_i |F_i|^2)`.

Validação: preserva os hashes T03, T04 e da regressão T05, confirma a relação sqrt(2) para as duas amplitudes RMS da varredura, remove ângulos espúrios da partícula central da cadeia e mantém o escopo em N=3, Lmax=1, sem Modelo D.
