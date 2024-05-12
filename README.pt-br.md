# Biblioteca de Bessel: uma biblioteca em C++ com rotinas para calcular funções de Bessel de argumentos reais ou complexos

<p align="right"><a href="README.md">Read in English</a></p>

## Recursos

### <nobr>*J*<sub>n</sub>(*z*)</nobr>
  - **Função:** `bessel::cyl_j(n,z)`.
  - **Descrição:** função cilíndrica de Bessel da primeira espécie para uma ordem inteira *n* e um argumento real ou compexo *z*.
  - **Opções:** a função aceita um terceiro argumento booleano opcional chamado *warnings*: `bessel::cyl_j(n,z,warnings)`. O valor padrão para *warnings* é `true`, defina-o como `false` para desativar os avisos. 
  - **Implementação para argumentos reais:** a rotina retorna a implementação da biblioteca padrão do C++.
  - **Implementação para argumentos complexos:** a rotina é baseada no cálculo de séries ascendentes, séries semiconvergentes de Stokes, recorrências diretas e recorrências retrógradas.
  A implementação segue a mesma ideia das rotinas apresentadas no Capítulo 5 da <nobr>Ref. [[1](#referências)]</nobr> para Fortran.

### <nobr>*J*<sub>n</sub>'(*z*)</nobr>
  - **Função:** `bessel::cyl_j_diff(n,z)`.
  - **Descrição:** Derivada de <nobr>*J*<sub>n</sub>(*z*)</nobr> com relação à *z*.
  - **Opções:** a função aceita um terceiro argumento booleano opcional chamado *warnings*: `bessel::cyl_j_diff(n,z,warnings)`. O valor padrão para *warnings* é `true`, defina-o como `false` para desativar os avisos. 
  - **Implementação para argumentos reais:** a rotina usa uma relação de recorrência junto com a implementação da biblioteca padrão do C++.
  - **Implementação para argumentos complexos:** a rotina retorna `-bessel::cyl_j(1,z)` se <nobr>*n* = 0</nobr> ou usa `bessel::cyl_j(n,z)` numa relação de recorrência caso contrário.
    
### Novos recursos
  - Infelizmente não há previsão para novos recursos.

## Como usar

A bibliotea está em estilo *header-only* (apenas cabeçalho), ou seja, não é necessário compilá-la separadamente, você só precisa incluir o arquivo *<a href="bessel-library.hpp">bessel-library.hpp</a>* no seu projeto.
Veja <a href="usage-example.cpp">usage-example.cpp</a> como exemplo de uso.

## Validação

As rotinas foram escritas para atingir resultados com precisão simples para <nobr>|Im(*z*)| ≤ 21</nobr> quando comparada com a função BesselJ[*n*,*z*] da <nobr>Ref. [[2](#referências)]</nobr>.

## Autoria

Os códigos e rotinas foram desenvolvidos e são atualizados por <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licença

Este projeto está protegido sob a licença <a href="LICENSE">MIT License</a>.

## Referências

[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.

[2] Wolfram Alpha, acessado em 01 Setembro 2023, <https://www.wolframalpha.com/>.