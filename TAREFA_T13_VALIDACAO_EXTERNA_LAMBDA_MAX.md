# T13 — validação externa de \(\Lambda_{\max}\)

## Objetivo e separação cronológica

A T13 valida externamente, sem reajuste, o preditor mecanístico M1 congelado
na T12.3. A amostra externa contém somente clusters planares com \(N=6\) e
\(N=10\), até então não usados para calibrar M1. O trabalho é dividido em uma
Fase A cega, publicada antes de qualquer resposta do Modelo E, e uma Fase B de
execução e análise. Nenhuma regra pode ser alterada entre as fases.

## Regime físico

São preservados \(a=1\), \(E_0=1\), \(ka=0.1\), \(f_0=0\), o plano nodal,
esferas idênticas e a força completa de interação do Modelo E como referência.
As forças A e D da T08 são apenas auditorias e diagnósticos. A força D não
substitui a referência E.

## Preditores congelados

O modelo confirmatório é

\[
\widehat\varepsilon_{M1}
=4.4964255121671126\,\Lambda_{\max}^{1.3883601043764593},
\qquad
\widehat\varepsilon_{M1,\mathrm{safe}}
=2.5699703122019222\,\widehat\varepsilon_{M1}.
\]

O comparador secundário, que não participa da decisão, é

\[
\widehat\varepsilon_{P3}
=14.73950709797405\,\rho_1^{1.4226504975598322},
\qquad
\widehat\varepsilon_{P3,\mathrm{safe}}
=2.0464420079866286\,\widehat\varepsilon_{P3}.
\]

M2 permanece excluído por `UNSTABLE_COLLINEARITY`. As classificações usam
estritamente `conservative_prediction < tolerance` e
`observed_error < tolerance`, para tolerâncias 0.01, 0.05 e 0.10.

## Seleção cega

O universo é o `split == holdout` da T08, limitado a \(N\in\{6,10\}\), às
famílias `linear`, `compact` e `irregular` e a referências D confirmadas.
Selecionam-se quatro casos distintos em cada um dos seis estratos, minimizando
a soma das distâncias logarítmicas aos alvos

\[
(0.0031111241226691642,\ 0.011108933664494051,\
0.025457132710914911,\ 0.065350897425260762).
\]

Empates até \(10^{-15}\) são resolvidos pela tupla lexicograficamente menor.
O checksum nominal é a lista `EXPECTED_CASE_IDS` em
`src/acoustic_ms/external_validation.py`. São 24 IDs, quatro por estrato; as
previsões conservadoras M1 classificam cegamente 6, 12 e 18 casos como seguros
em 1%, 5% e 10%.

## Artefatos imutáveis da Fase A

Antes de qualquer solve E, foram criados:

- `results/data/t13_holdout_manifest.csv`;
- `results/data/t13_frozen_predictions.csv`;
- `results/data/t13_frozen_protocol.csv`.

Os hashes SHA-256 iniciais são:

```text
25d79db59d9dd6d52c5674d0a64fe2fea351cf213a0cdcd92b45845a9ecc2b38  results/data/t13_holdout_manifest.csv
581a748dca2e5d161890284fca673ed20f2a4fbcbc0ff356d5d31db6ec8ac9c2  results/data/t13_frozen_predictions.csv
eb1878e3425ede7a2b599fd20f63550d2fdb23d177264a043d443694907dc650  results/data/t13_frozen_protocol.csv
```

Esses arquivos não contêm força, erro ou coluna de resposta E e tornam-se
imutáveis no commit `chore: preregister T13 external validation`.

## Campanha Modelo E

Cada caso é calculado em \(L_{\max}=2,3,\ldots,13\), sem parada antes de 5.
A parada antecipada exige confirmação por duas mudanças sucessivas aplicáveis,
menores ou iguais a \(10^{-5}\), em total, interação,
external–scattered e scattered–scattered. Se a interação não estiver
confirmada em 13, a execução prossegue até sua confirmação ou até 21. Nenhuma
força não confirmada é extrapolada.

Antes do solve E, geometria, \(\Lambda_{\max}\), \(\rho_1\), A e D são
auditados contra a T08. Cada ordem registra solver balanceado, finitude,
condicionamento menor que 10, erros de fechamento e decomposição menores que
\(10^{-12}\), simetria planar e dimensão modal.

## Elegibilidade, métricas e diagnósticos

Um caso participa das métricas somente se a interação E estiver confirmada,
os diagnósticos forem aprovados e o erro A–E for aplicável. São reportados os
24 casos, os elegíveis, cada valor de \(N\), família e nível-alvo. Para M1 e
P3 calculam-se RMSE e MAE logarítmicos, fatores multiplicativos, frações em
fator 2 e 1.5, Spearman e viés logarítmico, sem ajuste.

As auditorias conservadoras registram falsos seguros e falsos inseguros por
tolerância, globalmente e para cada \(N\). Diagnósticos mecanísticos incluem
erro D–E, amplitudes normalizadas dos canais, reforço ou cancelamento vetorial,
resíduos de M1 e extrapolação em relação à faixa de desenvolvimento. Eles não
alteram o gate.

## Gate congelado

A suficiência exige pelo menos 20 elegíveis no total, 10 para cada \(N\),
diagnósticos aprovados, manifesto íntegro e preservação da Fase A e dos 70
artefatos anteriores. Falha de suficiência produz
`INCONCLUSIVE_T13_INSUFFICIENT_MODEL_E_CONVERGENCE` e
`HOLD_T14_MODEL_E_CONVERGENCE`.

Com suficiência, M1 passa apenas com zero falsos seguros nas três tolerâncias;
antivacuidade 3, 6 e 9; ambos os valores de \(N\) representados em 5% e 10%;
RMSE global e por \(N\) não superior a \(\ln 2\); ao menos 80% global e 75%
por \(N\) em fator 2; Spearman global de pelo menos 0.90; e protocolo
imutável. O resultado é exatamente `PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX`
ou `FAIL_T13_EXTERNAL_VALIDATION_LAMBDA_MAX`. P3 não pode resgatar M1.

## Artefatos e verificações finais

A Fase B publica o CSV bruto de convergência atomicamente e deriva forças,
resumo de casos, previsões reveladas, métricas, auditoria de limiares, gate e
uma figura científica. O modo `--audit-existing` recalcula oito casos
estratificados; `--analyze-only` não chama o solver E e deve ser byte-idêntico
em duas execuções no mesmo ambiente.

A aprovação operacional requer suíte completa sem warnings, 24 IDs exatos,
nenhuma resposta usada na seleção, preservação byte a byte da Fase A e dos 70
artefatos anteriores, CSVs finitos e determinísticos, ausência de caracteres
de controle indevidos, auditoria estratificada e inspeção visual. Uma eventual
reprovação científica de M1 é resultado válido, não falha de software.

## Limitações

As conclusões restringem-se às famílias e parâmetros congelados, a
\(N\in\{6,10\}\), \(ka=0.1\), clusters planares e ao Modelo E disponível.
Não há reajuste, busca de novo preditor, T13.1 ou início automático da T14.
