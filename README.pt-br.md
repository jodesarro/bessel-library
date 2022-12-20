# Funções de Bessel

<p align="right"><a href="README.md">Read in English</a></p>

Funções de Bessel: uma biblioteca em C++ com rotinas para calcular funções de Bessel de argumentos reais ou complexos.

## Recursos

### <nobr>*J*<sub>0</sub>(*z*)</nobr>
  - Função: `bessel::cyl_j0(z)`.
  - Descrição: Função cilíndrica de Bessel de ordem zero para argumentos reais e compexos.
  - Implementação: Para argumento real, a rotina retorna a implementação da biblioteca padrão do C++.
  Para argumento complexo, a rotina calcula a série ascendente quando <nobr>|*z*| ≤ 12</nobr> e a série semiconvergente de Stokes caso contrário seguindo a mesma ideia apresentada no Capítulo 5 de [[1](#referências)] para Fortran.
  Vale ressaltar que para <nobr>|*z*| > 12</nobr>, `bessel::cyl_j0(z)` depende fortemente da precisão das funções `std::sin(z)` e `std::cos(z)`.

### Novos recursos
  - Infelizmente não há previsão para novos recursos.

## Como usar

A bibliotea está em estilo *header-only* (apenas cabeçalho), ou seja, não é necessário compilá-la separadamente, você só precisa incluir o arquivo *<a href="bessel-library.hpp">bessel-library.hpp</a>* no seu projeto.
Veja <a href="usage-example.cpp">usage-example.cpp</a> como exemplo de uso.

## Validação

A função `bessel::cyl_j0(z)` foi comparada com os valores tabulados do Capítulo 5.9 de [[1](#referências)] e com os valores resultantes da função BesselJ[0,*z*] de [[2](#referências)] dispostos na tabela a seguir.

|                                      | De [[1](#referências)]                                      | De [[2](#referências)]                                                     | De `bessel::cyl_j0(z)`                                                            |
|--------------------------------------|-------------------------------------------------------------|----------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
|<nobr>*J*<sub>0</sub>(1+*i*0)</nobr>  |+7.65197687×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |+7.651976865579665×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |+7.65197686557966**6**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(5+*i*0)</nobr>  |-1.77596771×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |-1.775967713143383×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |-1.77596771314338**4**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(10+*i*0)</nobr> |-2.45935764×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>         |-2.459357644513483×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                 |-2.459357644513**713**×10<sup>-1</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(25+*i*0)</nobr> |+9.62667833×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+9.626678327595811×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+9.62667832759581**2**×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(50+*i*0)</nobr> |+5.58123277×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+5.581232766925181×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+5.581232766925181×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                        |
|<nobr>*J*<sub>0</sub>(100+*i*0)</nobr>|+1.99858503×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>         |+1.998585030422312×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                 |+1.99858503042231**1**×10<sup>-2</sup><br/>+*i*0×10<sup>0</sup>                    |
|<nobr>*J*<sub>0</sub>(4+*i*2)</nobr>  |-1.3787022 ×10<sup>0</sup> <br/>+*i*3.9054236×10<sup>-1</sup>|-1.378702234392483×10<sup>0</sup> <br/>+*i*3.905423570667093×10<sup>-1</sup>|-1.37870223439248**4**×10<sup>0</sup><br/>+*i*3.90542357066709**4**×10<sup>-1</sup>|
|<nobr>*J*<sub>0</sub>(20+*i*10)</nobr>|+1.5460268 ×10<sup>+3</sup><br/>-*i*1.0391216×10<sup>+3</sup>|+1.546026837210333×10<sup>+3</sup><br/>-*i*1.039121575995158×10<sup>+3</sup>|+1.546026837210333×10<sup>+3</sup><br/>-*i*1.039121575995158×10<sup>+3</sup>       |

## Licença

Este projeto está protegido sob a licença <a href="LICENSE">MIT License</a> e tem [@jodesarro]( https://github.com/jodesarro ) como seu principal autor.

## Referências
[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.<br/>
[2] Wolfram Alpha, acessado em 09 abril 2022, <https://www.wolframalpha.com/>.
