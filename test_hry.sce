function [x,y,res]=testHry(A,xa,ya,N)
    sx = size(xa);
    sy = size(ya);
    
    xaA = xa'*A;
    Aya = A*ya;
    xaAya = xa'*A*ya;
    
    res = [%t %t];

    for i = 1:N
        [x,y] = vygenerujStrategie(sx,sy);
        res = porovnejStrategie(xaAya,xaA,Aya,x,y);

        if ~and(res)
            disp(i)
            break
        end
    end
    
    if res(1)
        x = xa;
    end
    if res(2)
        y = ya;
    end
    if and(res)
        disp('Optimalnejsi strategie nebyla nalezena.')
    else
        disp('Byla nalezena optimalnejsi strategie. ')
        disp(res)
    end
endfunction

function [x,y]=vygenerujStrategie(sx,sy)
    x = genStrategii(sx);
    y = genStrategii(sy);
endfunction

function x=genStrategii(s)
    x = rand(s(1),s(2));
    x = x./sum(x(:));
endfunction

function res=porovnejStrategie(xaAya,xaA,Aya,x,y)
    // x'.A.y*
    a = x'*Aya;
    // x*.A.y*
    b = xaAya;
    // x*.A.y'
    c = xaA*y;

    res = [(a<=b) (b<=c)];
endfunction
	
