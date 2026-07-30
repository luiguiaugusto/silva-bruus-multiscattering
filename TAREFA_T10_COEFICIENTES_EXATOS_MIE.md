# T10 — Coeficientes exatos de Mie de uma esfera fluida

## 1. Objetivo e base

A T10 parte do commit `729a2d5` e implementa a resposta parcial exata de uma
esfera fluida, compressível, homogênea e sem perdas. O objetivo é quantificar
o erro dos coeficientes líderes de Rayleigh no domínio atual do projeto. A
T-matrix nova permanece isolada: ela não é conectada ao Modelo D, não altera
as forças A--D e não constitui o Modelo E.

## 2. Convenções e derivação

Mantêm-se a dependência temporal (e^{-i\omega t}), os harmônicos esféricos
complexos do projeto e (h_\ell^{(1)}=j_\ell+i y_\ell). Com

\[
x=ka,\qquad
\widetilde\rho=\frac{\rho_p}{\rho_0},\qquad
\widetilde\kappa=\frac{\kappa_p}{\kappa_0},
\]

\[
y=x\sqrt{\widetilde\rho\widetilde\kappa},\qquad
\beta=\sqrt{\frac{\widetilde\kappa}{\widetilde\rho}},
\]

as continuidades de pressão e velocidade radial na superfície fornecem

\[
s_\ell^{\mathrm{Mie}}=-
\frac{
\beta j_\ell(x)j_\ell'(y)-j_\ell(y)j_\ell'(x)
}{
\beta h_\ell^{(1)}(x)j_\ell'(y)
-j_\ell(y){h_\ell^{(1)}}'(x)
}.
\]

As propriedades são obtidas dos contrastes de Silva--Bruus por

\[
\widetilde\kappa=1-f_0,\qquad
\widetilde\rho=\frac{2+f_1}{2(1-f_1)},\qquad
\frac{c_p}{c_0}=\frac{1}{\sqrt{\widetilde\rho\widetilde\kappa}},
\]

para (f_0<1) e (-2<f_1<1). O ponto exato (f_1=1) seleciona, sem
clipping nem densidade artificial, o limite rígido

\[
s_\ell^{\mathrm{rigid}}=-\frac{j_\ell'(x)}{{h_\ell^{(1)}}'(x)}.
\]

## 3. Implementação e estabilidade

`src/acoustic_ms/mie_scattering.py` separa a conversão de materiais, a esfera
fluida, o limite rígido e o wrapper por contrastes. As funções usam diretamente
`scipy.special.spherical_jn` e `spherical_yn`, inclusive as derivadas, e
retornam um vetor complexo para (ell=0,\ldots,L_{\max}). Escalares não reais,
booleanos, valores não finitos e domínios não físicos são rejeitados. A
correspondência material exata retorna zeros bit a bit. Fora do intervalo da
campanha, qualquer (ka>0) numericamente suportado é aceito, mas a escolha de
(L_{\max}) suficiente pertence ao usuário.

## 4. Validação independente

O teste independente monta, para cada ordem, o sistema (2\times2)

\[
\begin{pmatrix}
h_\ell(x)&-j_\ell(y)\\
h_\ell'(x)&-\beta j_\ell'(y)
\end{pmatrix}
\begin{pmatrix}s_\ell\\A_\ell\end{pmatrix}
=-
\begin{pmatrix}j_\ell(x)\\j_\ell'(x)\end{pmatrix}
\]

com `numpy.linalg.solve`. Para dois materiais, três valores de (ka) e
(ell=0,\ldots,5), a diferença absoluta máxima foi
(9.36\times10^{-20}); o maior erro relativo, no monopolo quase casado e
mal escalado, foi (3.17\times10^{-12}). O maior resíduo de condição de
contorno da campanha foi (2.22\times10^{-16}), e o máximo defeito de
unitariedade

\[
\left|\operatorname{Re}s_\ell+|s_\ell|^2\right|
\]

foi (1.32\times10^{-23}).

Os testes também cobrem a correspondência material, as potências
((ka)^{2\ell+1}), a correção relativa de ordem ((ka)^2), as fórmulas
publicadas para (ell=1,3,5), o limite rígido direto e por densidades
crescentes, unitariedade e validação de entradas. Na sequência
(\widetilde\rho=10^8,10^{10},10^{12}), o erro RMS relativo ao limite rígido
caiu de (1.63\times10^{-3}) para (2.05\times10^{-4}) e
(1.47\times10^{-5}).

## 5. Campanha Mie--Rayleigh

`scripts/analyze_t10_mie_rayleigh.py` usa (f_0=0),
(f_1\in\{0.1,0.4,0.8,1\}), 101 pontos logarítmicos em
(10^{-3}\le ka\le10^{-1}) e (ell=0,\ldots,5). O CSV principal contém
2.424 linhas e o resumo contém 24 linhas. Em (ka=0.1), os erros complexos
relativos do dipolo são:

| (f_1) | (\varepsilon_{s_1}) |
|---:|---:|
| 0.1 | 0.0018021200413863723 |
| 0.4 | 0.0012045872800514213 |
| 0.8 | 0.0004219360873243320 |
| 1.0, rígido | 0.0030243505842077237 |

As inclinações assintóticas observadas para todas as ordens positivas ficam
próximas de 2, como esperado para a primeira correção relativa. Como
(f_0=0), o termo Rayleigh de (ell=0) é zero, enquanto Mie conserva uma
correção dinâmica superior. Seu erro relativo é registrado com erro absoluto
e flag de aplicabilidade física; esse canal é inativo pela simetria do atual
problema nodal.

## 6. Artefatos e reprodutibilidade

Foram gerados:

```text
results/data/t10_mie_rayleigh_validation.csv
results/data/t10_mie_rayleigh_summary.csv
results/figures/t10_mie_rayleigh_error.png
```

Duas execuções consecutivas no ambiente Python 3.12.3, NumPy 2.5.1, SciPy
1.18.0 e Matplotlib 3.11.1 devem ser byte-idênticas. Os hashes oficiais desta
árvore de revisão são registrados em `docs/HANDOFF.md`.

## 7. Critérios de aceite e limitações

A tarefa exige todos os testes anteriores e novos sem warnings, oráculo e
condições de contorno aprovados, limites de Rayleigh e rígido recuperados,
unitariedade, artefatos determinísticos, figura inspecionada e preservação de
todos os resultados anteriores.

\[
\boxed{
\text{coeficientes exatos de Mie}
\ne
\text{força coletiva completa}
}
\]

A T10 valida somente a T-matrix diagonal de uma esfera ideal sem perdas. Ela
não valida independentemente Silva--Bruus ou o Modelo D, não inclui absorção,
viscosidade, elasticidade sólida, paredes, formas não esféricas nem termos
`scattered--scattered`. A integração ao solver global e uma formulação de
força completa pertencem à T11.
