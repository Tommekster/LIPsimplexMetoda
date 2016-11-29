function [sTab,B]=testSimplTab()
    sTab = [1 0 0 0 0.6 -6.4 4.8 0; 0 1 0 0 0.2 -1.8 0.6 0; 0 0 1 0 0.4 -1.6 0.2 0; 0 0 0 1 0 1 0 1; 0 0 0 0 -0.4 -0.4 1.8 0];
    B = [1 2 3 4];
endfunction

function c=sTab2c(sTab)
    c = sTab($,1:$-1);
endfunction

function A=sTab2A(sTab)
    A = sTab(1:$-1,1:$-1);
endfunction

function b=sTab2b(sTab)
    b = sTab(1:$-1,$);
endfunction

function a_s=sTab2As(sTab,s)
    a_s = sTab(1:$-1,s);
endfunction

function s=vyberSblandem(c)
    // s = min{j : cj < 0}
    ran = 1:size(c(:),1);
    s = min(ran(c<0));
endfunction

function r=vyberRblandem(As,b)
    // r t,z Br = min{Bk : bk/aks = min{bj/ajs : ajs > 0}, aks > 0}
    // {bj/ajs : ajs > 0}
    ran = 1:size(b,1);
    ran = ran(As>0);
    bjBYajs = b(ran)./As(ran);
    // r pro ktere plati br/ars = min
    mza_r = ran(bjBYajs == min(bjBYajs));
    r = min(mza_r);
endfunction

function [r,s]=blandovoPravidlo(sTab,B)
    c = sTab2c(sTab);
    b = sTab2b(sTab);
    
    // s = min{j : cj < 0}
    s = vyberSblandem(c);
    
    // r t,z Br = min{Bk : bk/aks = min{bj/ajs : ajs > 0}, aks > 0}
    r = vyberRblandem(sTab2As(sTab,s),b);
endfunction

// cp = ct-ctb.iAb.A
// Ap = iAb.A
function res=isNeomezena(As)
    res = max(As) <= 0;
endfunction

function A=eliminuj(A,r,s)
    // pripravim si matici radkovych uprav
    d = size(A,1);
    T = eye(d,d);
    T(:,r) = -A(:,s)./A(r,s);
    T(r,r) = 1/A(r,s);
    
    // nyni provedeme radkove upravy v jednom kroku
    A = T*A;

    // upravime sloupec na presne hodnoty, kvuli numericke stabilite
    A(:,s) = 0;
    A(r,s) = 1;
endfunction

function [sTab,B,optimalni] = simplexovaMetoda(sTab,B)
    optimalni = %f;
    neomezena = %f;
    omezeni = 15;
    
    while ~(optimalni | neomezena) | omezeni > 0
        omezeni = omezeni - 1;
        c = sTab2c(sTab);
        if min(c) >= 0
            optimalni = %t;
        else
            s = vyberSblandem(c);
            if isNeomezena(sTab2As(sTab,s))
                neomezena = %t;
            else
                // vyber r podle Blanda
                r = vyberRblandem(sTab2As(sTab,s),sTab2b(sTab));
                
                // eliminuj pivotem ars
                sTab = eliminuj(sTab,r,s);
                
                // doupresni bazi
                B(r) = s;
                
                // ukaz co mas v hlave
                disp(sTab);
                disp(B);
                //input('press a key...');
            end
        end
    end
    
    if optimalni 
        disp('xb je optimalni reseni')
    else
        disp('uloha je neomezena')
    end
endfunction

function [sTab,B]=simplexKrok(sTab,B)
    c = sTab2c(sTab);
    // vyber s,r podle Blanda
    s = vyberSblandem(c);
    r = vyberRblandem(sTab2As(sTab,s),sTab2b(sTab));
    
    // eliminuj pivotem ars
    sTab = eliminuj(sTab,r,s);
    
    // doupresni bazi
    B(r) = s;
    
    // ukaz co mas v hlave
    disp(sTab);
    disp(B);
endfunction

function xb=vytahniXb(sTab,B)
    xb = zeros(size(sTab,2)-1,1); // pripravime nulovy vektor
    xb(B) = sTab(1:$-1,$); // spravna mista vyplnime cilsy
endfunction

function h = hodnotaUceloveFce(sTab)
    h = -sTab($,$);
endfunction

function [sTab,B]=prvniFazeRozUlohy(A,b)
    [m,n] = size(A);
    
    B = (n+1):(n+m);
    
    sTab = zeros(m+1,n+m+1);
    sTab(1:m,1:n) = A;
    sTab(1:m,B) = eye(m,m);
    sTab(1:m,$) = b;
    // sTab($,1:n) = c; // c pouzijeme jine: c=(0,...,0,1,...,1)
    // sTab($,B) = 1; // takovy vektor je treba jeste eliminovat
    sTab($,1:n) = (-1.*ones(1,m))*A;
    sTab($,$) = -ones(1,m)*b; // -h=-(ctb.xbb)=-(ct.xb) 
endfunction

function [subTab,B]=upravBaziDoRozsahu1n(sTab,B)
    [d,e] = size(sTab);
    m = d-1;
    n = e-d;
    // Myslim, ze uz muzu cast tabulky zahodit
    subTab = sTab(:,[1:n,n+m+1]); 
    
    // opravit bazi
    ran = 1:(n+m);
    ranN = 1:n;
    spatnaBaze = ran(B>n);
    for r = spatnaBaze
        // najdu v radku jiny nenulovy prvek a tim zeliminuji
        //nenula = ranN(subTab(r,1:n) <> 0);
        // budu tim delit, tak vezmu ten nejvetsi; num. stab. doufam
        [v,s] = min(abs(subTab(r,1:n)));
        // konecne eliminuji pivotem
        subTab = eliminuj(subTab,r,s);
        // zaznamenam do baze
        B(r) = s; 
    end
endfunction

function [subTab,B] = rozsirenaUloha(A,b)
    [sTab,B] = prvniFazeRozUlohy(A,b);
    [sTab,B,optimalni] = simplexovaMetoda(sTab,B);
    
    minimum = sTab($,$); // zajima me, jestli minimum je nula
    if minimum == 0
        disp('Neexistuje pripustne reseni puvodni ulohy');
        subTab = 0;
        B = 0;
    else
        // Pripustne reseni existuje
        [subTab,B] = upravBaziDoRozsahu1n(sTab,B);
    end
endfunction

function sTab=tabulkaPuvodniUlohy(subTab,c,B)
    [d,e] = size(subTab);
    m = d-1;
    n = e-1;
    
    sTab = subTab;
    sTab($,1:n) = c;
    sTab($,$) = 0; // -h=-(ctb.xbb)=-(ct.xb) 
    //disp(sTab);
    for r=1:m
        //disp([r,B(r)]);
        sTab = eliminuj(sTab,r,B(r));
    end
endfunction

function [A,b,c]=testZadani()
    A = [-1,3,2,-1,-1,3,0,1,0,0;0,1,2,0,3,1,-1,0,0,0;
    -1,0,-1,2,1,2,2,-1,3,2;1,2,-1,2,-1,1,3,0,1,1;3,-1,0,3,0,3,1,0,1,0];
    b = [1;2;3;4;5];
    c = -[4,5,1,2,4,2,4,1,3,4];
endfunction
