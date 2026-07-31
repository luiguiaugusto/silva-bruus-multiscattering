# T14 — extrapolação de escala do critério congelado em \(\Lambda_{\max}\)

## Objetivo

A T14 valida externamente, sem recalibração, o preditor M1 congelado na T12.3
e aprovado na T13. A amostra contém \(N=15\) e \(N=28\), três famílias
planares e quatro níveis de acoplamento. O erro observado é sempre o erro RMS
normalizado entre o Modelo A e a força de interação completa do Modelo E.

## Separação cronológica

A fase A publica manifesto, previsões e protocolo antes de qualquer solução do
Modelo E nos 24 casos. Seus quatro CSVs são então imutáveis. Somente após o
commit `chore: preregister T14 lambda max scale-out` e seu push a fase B pode
executar a campanha. A conclusão científica é publicada em um segundo commit,
`feat: scale out frozen lambda max criterion`.

## Protocolo congelado

São fixados \(a=E_0=1\), \(ka=0.1\), \(f_0=0\) e \(f_1=0.8\). O modelo
confirmatório e sua margem conservadora são

\[
\widehat\varepsilon_{M1}=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\qquad
\widehat\varepsilon_{M1,\mathrm{safe}}
=2.5699703122019222\,\widehat\varepsilon_{M1}.
\]

P3 é apenas diagnóstico:

\[
\widehat\varepsilon_{P3}=14.73950709797405\,\rho_1^{1.4226504975598322},
\qquad
\widehat\varepsilon_{P3,\mathrm{safe}}
=2.0464420079866286\,\widehat\varepsilon_{P3}.
\]

A segurança usa comparações estritas nas tolerâncias 1%, 5% e 10%. M2 permanece
excluído por colinearidade instável. Nenhum desses coeficientes pode ser
alterado após a revelação.

## Amostra cega

As famílias linear, compacta triangular e irregular determinística são
normalizadas para distância mínima unitária, centralizadas e escaladas por

\[
d=\left(\frac{|f_1|S_{N,g}}{\Lambda_{\max}^{\mathrm{target}}}\right)^{1/3}.
\]

Os alvos são `0.0031111241226691642`, `0.011108933664494051`,
`0.025457132710914911` e `0.065350897425260762`. A irregular aplica uma
perturbação determinística de amplitude 0.15 baseada nas partes fracionárias de
\((i+1)\sqrt{2}\) e \((i+1)\sqrt{3}\), seguida de nova centralização e
normalização. A ordem nominal é tamanho, família `linear`, `compact`,
`irregular`, e nível 1–4. As contagens cegas M1 são exatamente 6, 12 e 18.

## Campanha Modelo E

Cada caso é resolvido nas ordens consecutivas \(L=2,\ldots,13\), com parada
antecipada a partir de \(L=5\) somente quando força total, interação,
external–scattered e scattered–scattered possuem confirmação em duas etapas.
O padrão é um trabalhador. O CSV bruto é publicado atomicamente e o cache
recuperável fica fora do repositório. Seis casos estratificados são
recalculados para auditoria.

## Elegibilidade e análise

Um caso é elegível quando a interação foi confirmada, os diagnósticos numéricos
passaram e o denominador físico do erro é aplicável. Casos inelegíveis não
entram nas métricas nem nas classificações. A análise não importa nem chama o
Modelo E e é reproduzida apenas do manifesto, protocolo, previsões e CSV bruto.

## Gates

A suficiência exige pelo menos 20 casos, 10 por tamanho, 6 por família e 5 por
nível, todos os casos cegamente seguros elegíveis, integridade dos artefatos e
nenhum solve acima de \(L=13\). Falha de suficiência retorna uma decisão
inconclusiva, nunca aprovação ou reprovação científica.

Com suficiência aprovada, M1 passa somente com zero falso seguro, contagens
6/12/18 preservadas, RMSE log global e por tamanho até \(\ln 2\), pelo menos
80% global e 75% por tamanho dentro de fator 2, Spearman global mínimo 0.90 e
protocolo imutável. P3 não participa da decisão.

## Integridade, determinismo e limitações

Todos os artefatos anteriores são congelados por caminho, tamanho e SHA-256.
A fase A e a análise da fase B devem ser byte-idênticas em duas execuções no
mesmo ambiente. O estudo é restrito a duas novas escalas, três formas fixas,
plano nodal, partículas idênticas e parâmetros canônicos. Um `PASS` valida
somente essa transferência; um `FAIL` é resultado científico publicável e não
autoriza recalibração nesta tarefa. T15 não é iniciada automaticamente.

## Resultado executado

A fase A foi publicada no commit
`6520173359b29cffa3a3d6432cefafcf17310f69` antes de qualquer solve E da
amostra. Os quatro artefatos cegos permaneceram byte-idênticos. A campanha
sequencial produziu 162 linhas de ordem, confirmou todos os quatro canais nos
24 casos e não excedeu L=11. Todos os casos foram elegíveis.

M1 obteve RMSE log global 0.343011370242051, fração dentro de fator 2 igual a
1.0 e Spearman 0.923685071658958. Houve zero falso seguro nas três tolerâncias
e as contagens 6/12/18 foram preservadas. A decisão literal foi
`PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX`, com próximo gate
`GO_T15_SYNTHESIS_AND_MANUSCRIPT`. Nenhum coeficiente foi recalibrado e T15 não
foi iniciada.
