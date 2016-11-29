exec('funkce.sce');

disp('Necht je zadana simplexova tabulka')
[A,b,c]=testZadani()

disp('Prubeh simplexove metody u pomocne ulohy')
[subTab,B] = rozsirenaUloha(A,b) 

disp('Simplexova tabulka pro puvodni ulohu')
sTab=tabulkaPuvodniUlohy(subTab,c,B)

disp('Prubeh simplexove metody pro puvodni ulohu')
[sTab,B] = simplexovaMetoda(sTab,B);

disp('Vysledek simplexove metody')
xb = vytahniXb(sTab,B)
h = hodnotaUceloveFce(sTab)
