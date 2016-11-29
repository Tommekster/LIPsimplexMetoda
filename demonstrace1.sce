exec('funkce.sce');

disp('Necht je zadana simplexova tabulka')
[sTab,B] = testSimplTab()

disp('Prubeh simplexove metody')
[sTab,B] = simplexovaMetoda(sTab,B);

disp('Vysledek simplexove metody')
xb = vytahniXb(sTab,B)
h = hodnotaUceloveFce(sTab)
