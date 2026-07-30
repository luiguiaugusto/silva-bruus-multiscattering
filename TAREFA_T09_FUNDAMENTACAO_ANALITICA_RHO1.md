# T09 — Fundamentação analítica de \(\rho_1\)

**Status:** concluída
**Base imutável:** T08, commit `48e88f8`
**Objetivo:** demonstrar por que o raio espectral do operador dipolar
balanceado mede a intensidade do reespalhamento coletivo e por que a primeira
correção à aproximação pairwise deve ser aproximadamente linear em \(\rho_1\).

## 1. Escopo

Esta fundamentação é exata dentro das seguintes hipóteses:

- esferas idênticas e fixas;
- centros no plano nodal \(z=0\);
- \(f_0=0\);
- campo acústico linear;
- coeficiente dipolar líder de Rayleigh;
- truncamento do campo espalhado em \(L=1\);
- observável de força external--scattered já validado na T04.

Ela não converte o Modelo D em uma T-matrix exata, não inclui a força
scattered--scattered e não demonstra que os limiares empíricos da T08 sejam
universais.

## 2. Sistema balanceado

O sistema físico de múltiplo espalhamento é

\[
(\mathbf I-\mathbf D\mathbf U)\mathbf s
=
\mathbf D\mathbf a_{\mathrm{ext}},
\]

onde \(\mathbf D\) contém os coeficientes de espalhamento de uma esfera e
\(\mathbf U\) translada o campo emitido por uma partícula para a base local de
outra. Definindo

\[
\mathbf s=\mathbf D^{1/2}\mathbf q,
\qquad
\mathbf b=\mathbf D^{1/2}\mathbf a_{\mathrm{ext}},
\]

obtemos

\[
(\mathbf I-\mathbf K_b)\mathbf q=\mathbf b,
\qquad
\boxed{\mathbf K_b=\mathbf D^{1/2}\mathbf U\mathbf D^{1/2}}.
\]

O balanceamento não introduz física nova. Ele apenas retira escalas
multipolares muito diferentes da matriz numérica. O operador \(\mathbf K_b\)
é a realimentação entre um evento de recepção e um novo evento de
espalhamento.

## 3. Redução exata no plano nodal em \(L=1\)

No plano nodal, a simetria de reflexão conserva os modos com
\(\ell+m\) ímpar. Para \(L=1\), \(f_0=0\) e o campo externo
\(\sin(kz)\), resta exatamente um canal por partícula:

\[
(\ell,m)=(1,0).
\]

O coeficiente dipolar líder é

\[
s_1=i\frac{f_1}{6}(ka)^3.
\]

Para dois centros separados por um vetor no plano \(xy\), o elemento de
translação desse canal é

\[
C(kr)
=
h_0^{(1)}(kr)+h_2^{(1)}(kr).
\]

Usando as formas fechadas dos Hankel esféricos,

\[
h_0^{(1)}(x)+h_2^{(1)}(x)
=
-3\frac{e^{ix}(x+i)}{x^3}.
\]

Como as esferas são idênticas, o produto das duas raízes principais de
\(s_1\) é \(s_1\). Logo, para \(i\ne j\),

\[
\begin{aligned}
(K_b)_{ij}
&=
s_1\left[
h_0^{(1)}(kr_{ij})+h_2^{(1)}(kr_{ij})
\right] \\
&=
\boxed{
\frac{f_1}{2}
\left(\frac{a}{r_{ij}}\right)^3
e^{ikr_{ij}}(1-ikr_{ij})
},
\end{aligned}
\]

e

\[
(K_b)_{ii}=0.
\]

Portanto, o operador \(L=1\) usado para calcular \(\rho_1\) pode ser
representado exatamente por uma matriz complexa \(N\times N\), sem montar a
base multipolar completa.

## 4. Origem da escala inverso-cubo

Para \(x=kr\ll1\),

\[
e^{ix}(1-ix)
=
1+\frac{x^2}{2}
+i\frac{x^3}{3}
-\frac{x^4}{8}
+O(x^5).
\]

Não existe correção linear em \(kr\). Assim,

\[
(K_b)_{ij}
=
\frac{f_1}{2}
\left(\frac{a}{r_{ij}}\right)^3
\left[
1+O((kr_{ij})^2)
\right].
\]

Esse resultado fornece a origem analítica da escala usada anteriormente:

\[
\eta=|f_1|\left(\frac{a}{d_{\min}}\right)^3.
\]

No dímero,

\[
\boxed{
\rho_1
=
\frac{|f_1|}{2}
\left(\frac{a}{d}\right)^3
\sqrt{1+(kd)^2}
}
\]

e, no limite próximo,

\[
\rho_1\longrightarrow\frac{\eta}{2}.
\]

Para um cluster geral,

\[
\rho_1\leq\|\mathbf K_b\|_\infty
=
\frac{|f_1|}{2}
\max_i\sum_{j\ne i}
\left(\frac{a}{r_{ij}}\right)^3
\sqrt{1+(kr_{ij})^2}.
\]

No limite quase estático, a matriz torna-se real e simétrica. Para
\(f_1>0\), ela também é não negativa, de modo que o teorema de
Perron--Frobenius fornece

\[
\min_i R_i
\leq
\rho_1^{\mathrm{nf}}
\leq
\max_i R_i
=
\frac{\Lambda_{\max}}{2},
\]

com

\[
R_i=\frac{|f_1|}{2}
\sum_{j\ne i}\left(\frac{a}{r_{ij}}\right)^3.
\]

Logo, \(\eta\) mede somente a ligação mais próxima,
\(\Lambda_{\max}\) fornece um limite coletivo por soma de linhas e
\(\rho_1\) mede os modos coletivos efetivos do grafo completo.

## 5. Série de Neumann e significado de reespalhamento

O sistema balanceado também pode ser escrito como

\[
\mathbf q=\mathbf b+\mathbf K_b\mathbf q.
\]

A solução formal é

\[
\boxed{
\mathbf q
=
\sum_{p=0}^{\infty}\mathbf K_b^p\mathbf b
}
\]

A série matricial de Neumann que sustenta essa expressão converge, para toda
fonte \(\mathbf b\), se e somente se

\[
\boxed{\rho(\mathbf K_b)<1}.
\]

Cada potência possui interpretação física:

- \(p=0\): cada esfera espalha apenas o campo externo;
- \(p=1\): um reespalhamento adicional;
- \(p=2\): dois reespalhamentos adicionais;
- e assim sucessivamente.

Após truncar em \(P\),

\[
\mathbf q-\mathbf q^{(P)}
=
\mathbf K_b^{P+1}
(\mathbf I-\mathbf K_b)^{-1}\mathbf b.
\]

Se uma norma subordinada satisfizer \(\|\mathbf K_b\|<1\),

\[
\|\mathbf q-\mathbf q^{(P)}\|
\leq
\frac{\|\mathbf K_b\|^{P+1}}
{1-\|\mathbf K_b\|}
\|\mathbf b\|.
\]

Para um operador diagonalizável,

\[
\|\mathbf K_b^p\|
\lesssim
\kappa(\mathbf V)\rho_1^p,
\]

onde \(\mathbf V\) contém os autovetores. Portanto, \(\rho_1\) controla a
taxa assintótica, enquanto a não normalidade pode modificar a amplificação
em ordens finitas.

## 6. Por que o erro pairwise começa em \(O(\rho_1)\)

No Modelo \(L=1\), a força external--scattered é uma transformação real
linear dos coeficientes espalhados. Denotando-a por \(\mathcal R\),

\[
\mathbf F^{D_1}
=
\mathcal R\!\left[
\mathbf D^{1/2}
(\mathbf I-\mathbf K_b)^{-1}\mathbf b
\right].
\]

O termo de ordem zero reproduz a soma pairwise de Silva--Bruus:

\[
\mathbf F^A
=
\mathcal R\!\left[
\mathbf D^{1/2}\mathbf b
\right].
\]

Consequentemente,

\[
\mathbf F^{D_1}-\mathbf F^A
=
\mathcal R\!\left[
\mathbf D^{1/2}
\sum_{p=1}^{\infty}\mathbf K_b^p\mathbf b
\right].
\]

O primeiro termo omitido é linear em \(\mathbf K_b\). Na ausência de
cancelamentos excepcionais da força de referência,

\[
\boxed{\varepsilon_A=O(\rho_1)}.
\]

Isso explica por que o ajuste empírico da T08 produziu

\[
\varepsilon_A
\simeq
2.635\,\rho_1^{1.109},
\]

com expoente próximo de um. A álgebra explica a ordem dominante; ela não
determina o prefator \(2.635\), pois esse número depende da geometria, da
projeção da força, das fases e da normalização pelo valor de referência.

A T08 compara \(A\) com \(D_L\), e não somente com \(D_1\). Assim, a relação
com a força multipolar convergida permanece uma extensão empírica apoiada na
dominância dipolar em \(ka=0.1\), não um limite rigoroso para todos os
multipolos.

## 7. Hierarquia de muitos corpos

O observável de força já contém uma propagação terminal entre a esfera que
espalha e a esfera que recebe a força. Cada fator adicional de
\(\mathbf K_b\) acrescenta outra aresta ao caminho de reespalhamento.

Um caminho conectado que visita \(n\) partículas distintas precisa de pelo
menos \(n-1\) arestas. Como uma delas já está no observável terminal, a
primeira contribuição irredutível de \(n\) corpos exige

\[
p_{\min}=n-2
\]

fatores de \(\mathbf K_b\). Logo,

\[
\boxed{
\Phi^{(3)}=O(\rho_1),
\qquad
\Phi^{(4)}=O(\rho_1^2)
}.
\]

Essa contagem explica as leis observadas na T06.1:

\[
Y_3\sim\Lambda,
\qquad
Y_4\sim\Lambda^2.
\]

O baseline pairwise corrigido \(B_1\) reabsorve todos os passeios confinados
a um dímero. A diferença \(D_1-B_1\) retém os passeios que visitam ao menos
três partículas e também começa, em geral, em primeira ordem. Isso é
consistente com o ajuste

\[
\varepsilon_B\simeq0.651\,\rho_1^{1.014}.
\]

## 8. Por que \(\rho_1\) contém mais informação

\(\eta\) conhece apenas \(d_{\min}\). Dois clusters com a mesma menor
separação podem ter conectividades inteiramente diferentes.

\(\Lambda_{\max}\) soma magnitudes inverso-cubo incidentes sobre a partícula
mais acoplada, mas não resolve:

- fases de propagação;
- sinais dos modos coletivos;
- ciclos fechados;
- interferência entre caminhos;
- distribuição global dos autovetores.

\(\rho_1\) é calculado do operador completo e identifica a maior taxa de
realimentação entre os modos coletivos permitidos. Isso justifica sua
superioridade modesta, mas sistemática, na validação cruzada da T08.

## 9. Não normalidade e limites do raio espectral

\(\mathbf K_b\) é complexo e simétrico, mas não é necessariamente
Hermitiano nem normal. Portanto,

\[
\rho_1\leq\|\mathbf K_b\|_2
\]

e o raio espectral, sozinho, não limita toda amplificação transitória.

Isso não invalida seu papel:

- \(\rho_1<1\) é a condição necessária e suficiente para convergência da
  série de Neumann em dimensão finita;
- \(\rho_1\) fornece a taxa assintótica;
- normas e medidas de não normalidade complementam o diagnóstico quando há
  autovetores mal condicionados.

Nos 312 casos da T08,

\[
\max\rho_1=0.2544601332,
\]

\[
\max\frac{\|\mathbf K_b\|_2}{\rho_1}
=
1.0174745603,
\]

\[
\max\frac{\|\mathbf K_b\|_\infty}{\rho_1}
=
1.4935853869.
\]

Assim, todos os casos estão no domínio convergente e o operador mostrou-se
próximo de normal na medida mais diretamente relevante para amplificação
espectral. A medida normalizada do comutador atingiu \(0.1602877434\), mas
isso ocorreu em acoplamento fraco e não produziu separação importante entre
\(\|\mathbf K_b\|_2\) e \(\rho_1\).

## 10. Auditoria reproduzível

A implementação independente da expressão fechada foi comparada com o
\(\rho_1\) congelado da T08 em todas as 312 configurações:

\[
\max
\left|
\rho_1^{\mathrm{analítico}}
-\rho_1^{\mathrm{T08}}
\right|
=
4.1633363423\times10^{-16}.
\]

As identidades de Hankel e a expansão de campo próximo foram verificadas
simbolicamente com SymPy. Não foi necessário usar Wolfram Mathematica.

No subconjunto com \(d_{\min}/a=2.1\), a diferença relativa máxima entre o
raio espectral exato do modelo dipolar e sua aproximação quase estática foi
\(3.985\%\). Em todo o domínio ela alcançou \(35.915\%\), sobretudo nas
configurações mais diluídas, em que \(kr\) já não é pequeno, mas o
acoplamento absoluto é muito fraco. Por isso, a lei inverso-cubo fornece a
escala física, enquanto o \(\rho_1\) usado no critério deve continuar sendo
calculado com o fator retardado completo.

A série de Neumann foi comparada com a solução direta em três casos
sentinela:

- \(N=4\), \(\rho_1=0.00599\);
- \(N=6\), \(\rho_1=0.0672\);
- \(N=10\), \(\rho_1=0.254\).

O erro das somas parciais decaiu monotonicamente até o piso de
arredondamento e acompanhou a taxa \(\rho_1^P\).

## 11. O que a T09 demonstra

Fica demonstrado que:

1. \(\rho_1\) é o raio espectral de um operador físico de reespalhamento, não
   uma métrica escolhida apenas por regressão;
2. a escala \(|f_1|(a/r)^3\) decorre analiticamente do canal dipolar;
3. a primeira correção pairwise é de primeira ordem no acoplamento;
4. contribuições conectadas de três e quatro corpos começam,
   respectivamente, nas ordens um e dois;
5. a série coletiva é convergente em todos os casos da T08;
6. a não normalidade é uma limitação conceitual real, mas pequena no domínio
   amostrado.

Não fica demonstrado que:

1. o prefator \(2.635\) ou o expoente \(1.109\) sejam universais;
2. \(\rho_1\) forneça um limite rigoroso do erro da força multipolar completa;
3. os limiares de 1%, 5% e 10% se transfiram para outros \(ka\), contrastes,
   geometrias ou modelos de força;
4. \(\rho_1<1\) seja suficiente para tornar Silva--Bruus preciso — essa
   condição garante convergência do reespalhamento, não pequenez do erro.

## 12. Artefatos

```text
src/acoustic_ms/rho_foundation.py
scripts/analyze_t09_rho_foundation.py
tests/test_t09_rho_foundation.py
tests/test_t09_artifacts.py
results/data/t09_operator_audit.csv
results/data/t09_neumann_convergence.csv
results/data/t09_analytic_summary.csv
results/figures/t09_rho_foundation.png
```

Execução:

```bash
.venv/bin/python scripts/analyze_t09_rho_foundation.py
.venv/bin/python -m pytest -q -W error
```

Nenhum CSV ou resultado das T01--T08 é reescrito.
