# T07 — Modelo D e convergência multipolar

## Objetivo e escopo

Implementar o Modelo D como multiple scattering global no plano nodal com
ordem multipolar variável, sem alterar os Modelos A–C nem a definição
external–scattered da força. O estudo cobre dímeros e as seis geometrias
canônicas com (N\leq4), (ka=0.1), esferas idênticas e sem viscosidade,
streaming ou dinâmica.

## Espalhamento e base

O monopolo permanece

\[
s_0=-\frac{i}{3}f_0(ka)^3,
\]

e, para \(\ell\geq1\), usa-se o termo Rayleigh dominante

\[
s_\ell=i\frac{3\ell f_1}{(2\ell-1)!!(2\ell+1)!!
[2(2\ell+1)-(\ell-1)f_1]}(ka)^{2\ell+1}.
\]

Esses coeficientes não são uma T-matrix exata em frequência arbitrária. A
ordenação pública contém os \((L_{\max}+1)^2\) modos. A base planar ativa
retém \(\ell+m\) ímpar e remove coeficientes de espalhamento exatamente
nulos; modos inativos são reconstruídos como zeros.

## Solver balanceado e força

A equação física e sua forma balanceada são

\[
(\mathbf I-\mathbf D\mathbf U)\mathbf s=\mathbf D\mathbf a_{\rm ext},
\]

\[
(\mathbf I-\mathbf D^{1/2}\mathbf U\mathbf D^{1/2})\mathbf q
=\mathbf D^{1/2}\mathbf a_{\rm ext},\qquad
\mathbf s=\mathbf D^{1/2}\mathbf q.
\]

O resíduo é sempre avaliado na primeira equação. A força reexpande todas as
fontes até alvo \(\ell=2\), exclui o campo próprio e reutiliza exatamente as
combinações de \(b_{2,-1}\) e \(b_{2,1}\) da T04. Ordem multipolar e número de
reespalhamentos são conceitos distintos: o sistema linear já ressoma os
reespalhamentos admitidos pela base.

## Validação

O Modelo D em \(L_{\max}=1\) deve reproduzir o Modelo C. Um sistema reduzido
independente do dímero valida \(L_{\max}=1,3,5\). A Eq. (30) é comparada ao
ramo estritamente reduzido usado em sua derivação nos valores
\(ka=0.1,0.05,0.025\), verificando a tendência assintótica sem tratar a fórmula
fechada como oráculo do solver planar geral.

O mapa do dímero usa \(d/a\in\{2,2.05,2.1,2.5,3\}\),
\(f_1\in\{0.1,0.4,0.8,1\}\) e \(L_{\max}\in\{1,3,5,7,9\}\), com \(L=11\)
somente para confirmação direcionada. A convergência exige duas diferenças
sucessivas abaixo de \(10^{-3}\).

Nas seis geometrias canônicas, todos os subconjuntos são resolvidos na mesma
ordem e os termos conectados são obtidos por inclusão–exclusão vetorial. As
cadeias podem exigir \(L=11\), e somente o quarteto linear pode avançar a
\(L=13\).

## Artefatos e aceite

São gerados três CSVs e duas figuras `t07_*`. A tarefa é aceita quando os
testes anteriores e novos passam sem warnings, as execuções consecutivas são
byte-idênticas, os artefatos anteriores permanecem inalterados, as figuras
são legíveis e nenhum arquivo científico protegido é modificado.

## Limitações

O estudo permanece no regime de Rayleigh, plano nodal, força
external–scattered, \(N\leq4\), esferas idênticas e fixas. Não inclui
T-matrix exata, multipolaridade de Mie, scattered–scattered, viscosidade,
streaming, paredes, torque ou dinâmica.
