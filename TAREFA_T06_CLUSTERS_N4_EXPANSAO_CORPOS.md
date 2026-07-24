# T06 — clusters com \(N=4\) e expansão conectada até quatro corpos

## Objetivo e regime

A T06 estende a comparação A/B/C da T05 a quartetos planares, preservando o solver Rayleigh com \(L_{\max}^{\mathrm{scatter}}=1\) e avaliação local até \(\ell=2\). O objetivo é separar a correção coletiva em contribuições irredutíveis de três corpos já definidas nos trímeros e uma contribuição genuinamente irredutível de quatro corpos. Modelo D, multipolos superiores, dinâmica, torque e termos espalhado--espalhado permanecem fora do escopo.

## Modelos e expansão por subconjuntos

O Modelo A soma forças pairwise Silva--Bruus. O Modelo B soma os seis pares isolados resolvidos pelo mesmo Modelo C Rayleigh. O Modelo C resolve simultaneamente as quatro partículas, com 16 coeficientes complexos. Para cada trímero \(T\) que contém a partícula \(i\),

\[
\boldsymbol{\Phi}_i^{(3)}(T)=\mathbf F_i^C(T)-\sum_{j\in T,\,j\ne i}\mathbf F_i^C(\{i,j\}).
\]

Os trímeros seguem a ordem `(0, 1, 2)`, `(0, 1, 3)`, `(0, 2, 3)`, `(1, 2, 3)`. Define-se

\[
\mathbf F_i^{(\le3)}=\mathbf F_i^B+\boldsymbol{\Phi}_{i,\Sigma}^{(3)},
\qquad
\boldsymbol{\Phi}_i^{(4)}=\mathbf F_i^C(Q)-\mathbf F_i^{(\le3)}.
\]

A forma fechada auditada é

\[
\boldsymbol{\Phi}_i^{(4)}=\mathbf F_i^C(Q)
-\sum_{T\ni i,\,|T|=3}\mathbf F_i^C(T)
+\sum_{j\ne i}\mathbf F_i^C(\{i,j\}).
\]

Todas as identidades são vetoriais e usam componentes cartesianas assinadas. Para \(N=4\), \(C-B\) contém tanto \(\boldsymbol{\Phi}_{\Sigma}^{(3)}\) quanto \(\boldsymbol{\Phi}^{(4)}\).

## API e geometrias

`decompose_nodal_quartet` aceita exatamente uma matriz finita `(4, 3)`, exige `lmax=1`, reutiliza `compare_nodal_force_models`, preserva a comparação completa e as quatro comparações de trímeros e retorna os campos A/B/C, correções de dois corpos, correção coletiva, contribuições por trímero, soma de três corpos, reconstrução até três corpos e termo de quatro corpos. As geometrias canônicas são a cadeia linear, o quadrado e o quadrilátero irregular fixo, todos planares, centrados e parametrizados por \(d_{\min}\).

## Testes, produção e aceite

O oráculo escalar de teste usa diretamente funções esféricas da SciPy e `numpy.linalg.solve`; valida pares, trímeros, quarteto, coeficientes \(s_{10}\), modos proibidos e termos conectados. A cobertura inclui entradas, incorporação lexicográfica, identidades recursiva e fechada, permutação, translação, rotação, escalamento, energia, contrastes, cadeia, quadrado, irregular e acoplamento fraco.

O validador usa \(a=E_0=1\), \(ka=0.1\), \(f_0=0\), quatro valores de \(f_1\) e 160 distâncias em \([2.1,10]\), totalizando 1.920 configurações. Produz CSV de 12 regressões, CSV da varredura e duas figuras. O aceite exige testes sem warnings, regressões canônicas, hashes determinísticos, artefatos anteriores intactos, ausência de caracteres de controle e nenhum arquivo científico protegido alterado.
