# LIP: simplexová metoda

> Zápočtový program do předmětu Lineární programování (LIP)

## Použití

Algoritmus je napsán pro prostředí ***Scilab***, které je velice podobné 
prostředím *GNU Octave* nebo *MATLAB*.

To hlavní z lineárního programování je obsaženo v jednotlivých funkcích 
souboru *funkce.sce*. Funkčnost a samotné použítí lze vyčíst ze skriptů:

* demostrace1.sce
* demostrace2.sce
* hra.sce

Pro spuštění zmíněných skriptů stačí:

1. otevřít *Scilab* v adresáři, který obsahuje tyto soubory
2. spustit příkaz `exec('demonstrace1.sce')` (bez středníku na konci ... je vidět průběh)

### Spuštění *Scilab* na FJFI

Nenalezl jsem *Scilab* na *kmlinux*, ale je možné se připojit na
*newton.fjfi.cvut.cz*, který běží na KFE.
Jinak je možnost si *Scilab* [stáhnout][scilab].

## Aplikace algoritmu na maticové hry
Tato simplexová metoda v základu řeší úlohu `min{c.x : Ax=b, x>=0}`,
avšak vhodnou formulací lze na tuto úlohu převést i úlohy následující:
``min{c.x : Ax>=b, x>=0}`` a
``max{c.x : Ax<=b, x>=0}``.
Způsob jak tyto úlohy přeformulovat lze vykoukat ze skriptu *hra.sce*, 
nebo také přímo od "zdroje" na prvních stránkách [wikiskript][wikiskripta].

Soubor *hra.sce* je ještě doplněn souborem *test_hry.sce*, který obsahuje 
funkce pro generování náhodných strategií a jejich porovnání. Slouží k
jednoduchému testu, jestli zadané strategie jsou optimální, nebo jestli 
náhodně vybere nějakou optimálnější. 
**Test se na konci skriptu *hra.sce* spustí.**

* * *

## Reference
1. J. Rohn: Lineární a nelineární programování (Učební text)
2. [Wikiskripta: 01LIP][wikiskripta]
3. [www.scilab.org][scilab]

[wikiskripta]: http://wikiskripta.fjfi.cvut.cz/wiki/index.php/01LIP "Wikiskripta: 01LIP"
[scilab]: http://www.scilab.org/ "Home - Scilab"
