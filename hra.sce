exec('funkce.sce');

// Zadani matice
A=[-1,3,-1,-2,1;0,-1,2,-1,1;1,1,3,2,-1;1,-2,2,-1,3;0,0,-1,3,-1]
//A = [0,1,-1;-1,0,1;1,-1,0]

// vyber parametru alpha: > 0 ^ > min(-a_ij)
alpha = max(0,max(-A(:))) + 0.1

// utvor matici E a vektor e
E = ones(A);
e = ones(size(A,2),1);

// matice s pruham, transponovana
Apt = (A + alpha*E)'

// chci nerovnost, takze delam AI.xp=b <=> Apt.x<=b
// tzn. zavedu mekke promenne
m = size(Apt,1);
AI = [Apt -eye(m,m)] // nerovnost: vetsi,rovno
c = [e' zeros(1,m)]; // mekke promenne nejsou v ucelove fci

// priprav simplexovou tabulku
[subTab,B] = rozsirenaUloha(AI,e)
sTab = tabulkaPuvodniUlohy(subTab,c,B)
//[m,n] = size(Apt);
//sTab = zeros(m+1,n+m+1);
//sTab(1:m,1:n) = -Apt;
//sTab($,1:n) = e'; 
//sTab(1:m,(n+1):(n+m)) = eye(m,m);
//B = (n+1):(n+m);

// simplexova metoda pustena na maticovou hru
[sTab,B] = simplexovaMetoda(sTab,B);
xsb = vytahniXb(sTab);
x0 = xsb(1:size(Apt,2)) // ozn. podle V22
xa = x0./(e'*x0) // x^asterisk

// pro jistotu prubnu i y
Ap = Apt';
m = size(Ap,1);
e = ones(size(A,1),1);
c = -[e' zeros(1,m)]; // hledam maximum
AI = [Ap eye(m,m)]; // nerovnost: mensi,rovno

[subTab,B] = rozsirenaUloha(AI,e)
sTab = tabulkaPuvodniUlohy(subTab,c,B)

[sTab,B] = simplexovaMetoda(sTab,B);
ysb = vytahniXb(sTab);
y0 = ysb(1:size(Ap,2)) // ozn. podle V22
ya = y0./(e'*y0) // y^asterisk

// zaver
disp('Optimalni strategie a cena hry')
xa,ya
omega = xa'*A*ya