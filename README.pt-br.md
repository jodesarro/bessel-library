# Biblioteca de Bessel: uma biblioteca em C++ com rotinas para calcular funções de Bessel de argumentos reais ou complexos

<p align="right"><a href="README.md">Read in English</a></p>

## Recursos

### <nobr>*J*<sub>n</sub>(*z*)</nobr>
  - **Função:** `bessel::cyl_j(n,z)`.
  - **Descrição:** função cilíndrica de Bessel da primeira espécie para uma ordem inteira *n* e um argumento real ou compexo *z*.
  - **Implementação para argumentos reais:** a rotina retorna a implementação da biblioteca padrão do C++.
  - **Implementação para argumentos complexos:** para <nobr>|*n*| ≤ 1</nobr>, a rotina calcula a série ascendente se <nobr>|*z*| ≤ 12</nobr> ou a série semiconvergente de Stokes caso contrário.
  Para <nobr>|*n*| > 1</nobr>, a rotina faz uso de uma recorrência progressiva se <nobr>*n* < |*z*|/4</nobr> ou uma recorrência regressiva caso contrário.
  A implementação segue a mesma ideia das rotinas apresentadas no Capítulo 5 da <nobr>Ref. [[1](#referências)]</nobr> para Fortran.

### <nobr>*J*<sub>n</sub>'(*z*)</nobr>
  - **Função:** `bessel::cyl_j_diff(n,z)`.
  - **Descrição:** Derivada de <nobr>*J*<sub>n</sub>(*z*)</nobr> com relação à *z*.
  - **Implementação para argumentos reais:** a rotina usa uma relação de recorrência junto com a implementação da biblioteca padrão do C++.
  - **Implementação para argumentos complexos:** a rotina retorna `-bessel::cyl_j(1,z)` se <nobr>*n* = 0</nobr> ou usa `bessel::cyl_j(n,z)` numa relação de recorrência caso contrário.
    
### Novos recursos
  - Infelizmente não há previsão para novos recursos.

## Como usar

A bibliotea está em estilo *header-only* (apenas cabeçalho), ou seja, não é necessário compilá-la separadamente, você só precisa incluir o arquivo *<a href="bessel-library.hpp">bessel-library.hpp</a>* no seu projeto.
Veja <a href="usage-example.cpp">usage-example.cpp</a> como exemplo de uso.

## Validação

As rotinas foram escritas para atingir resultados com precisão simples para <nobr>Im(*z*) ≤ 22</nobr>.

Para determinados parâmetros, as funções foram comparadas com a função BesselJ[*n*,*z*] da <nobr>Ref. [[2](#referências)]</nobr>. Os resultados estão dispostos na tabela a seguir.

|                                            | Da Ref. [[2](#references)]                                                   | Desta biblioteca                                                                    |
|--------------------------------------------|--------------------------------------------------------------------------------|--------------------------------------------------------------------------------------|
|<nobr>*J*<sub>0</sub>(3.2-*i*1.7)</nobr>    |-1.030282834527403×10<sup>0</sup><br/>  +*i*5.419144653871709×10<sup>-1</sup>   |-1.03028283452740**4**×10<sup>0</sup><br/>  +*i*5.41914465387170**5**×10<sup>-1</sup> |
|<nobr>*J*<sub>0</sub>(-15.1-*i*83.9)</nobr> |-9.079800270662861×10<sup>+34</sup><br/>-*i*7.605394317402608×10<sup>+34</sup>  |-9.079800270662**912**×10<sup>+34</sup><br/>-*i*7.6053943174026**56**×10<sup>+34</sup>|
|<nobr>*J*<sub>-1</sub>(-2.8+*i*7.7)</nobr>  |+1.405657382170826×10<sup>+2</sup><br/> +*i*2.583924513755478×10<sup>+2</sup>   |+1.40565738217082**8**×10<sup>+2</sup><br/> +*i*2.58392451375547**9**×10<sup>+2</sup> |
|<nobr>*J*<sub>3</sub>(4+*i*29)</nobr>       |+1.807884960557100×10<sup>+11</sup><br/>+*i*1.717728203190626×10<sup>+11</sup>  |+1.807884960557100×10<sup>+11</sup><br/>    +*i*1.717728203190626×10<sup>+11</sup>    |
|<nobr>*J*<sub>-14</sub>(42-*i*8)</nobr>     |-1.109381836657852×10<sup>+2</sup><br/> +*i*4.649254392531784×10<sup>+1</sup>   |-1.109381836657852×10<sup>+2</sup><br/>     +*i*4.64925439253178**5**×10<sup>+1</sup> |
|<nobr>*J*<sub>100</sub>(71+*i*30)</nobr>    |-1.873523256314857×10<sup>-4</sup><br/> +*i*2.688597190207697×10<sup>-6</sup>   |-1.8735232563148**47**×10<sup>-4</sup><br/> +*i*2.68859719020**8606**×10<sup>-6</sup> |
|<nobr>*J*<sub>0</sub>'(-11+*i*18)</nobr>    |-5.427270167814426×10<sup>+6</sup><br/> +*i*1.444852021146701×10<sup>+6</sup>   |-5.427270167814426×10<sup>+6</sup><br/>     +*i*1.44485202114670**2**×10<sup>+6</sup> |
|<nobr>*J*<sub>-9</sub>'(-6.6-*i*3.6)</nobr> |+1.608080761747265×10<sup>-1</sup><br/> -*i*1.505778965628776×10<sup>-1</sup>   |+1.60808076174726**4**×10<sup>-1</sup><br/> -*i*1.50577896562877**1**×10<sup>-1</sup> |
|<nobr>*J*<sub>47</sub>'(2.71+*i*3.14)</nobr>|-4.836547831119548×10<sup>-45</sup><br/>+*i*3.403725372134759×10<sup>-44</sup>  |-4.8365478311195**86**×10<sup>-45</sup><br/>+*i*3.40372537213475**6**×10<sup>-44</sup>|

## Autoria

Os códigos e rotinas foram desenvolvidos e são atualizados por <a href="https://www.researchgate.net/profile/Jhonas-de-Sarro">Jhonas O. de Sarro</a> ([@jodesarro]( https://github.com/jodesarro )).

## Licença

Este projeto está protegido sob a licença <a href="LICENSE">MIT License</a>.

## Referências

[1] S. Zhang and J. Jin, "Computation of Special Functions", Wiley, 1996.

[2] Wolfram Alpha, acessado em 01 Setembro 2023, <https://www.wolframalpha.com/>.