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
    bjBYajs = b(As>0)./As(As>0);
    // r pro ktere plati br/ars = min
    ran = 1:size(b,1);
    ran = ran(As>0);
    mza_r = ran(bjBYajs == min(bjBYajs));
    r = min(B(mza_r));
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
