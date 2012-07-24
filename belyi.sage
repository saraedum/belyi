def eindToFront(S):
    #   Die Eingabe sei eine absteigend sortierte Liste

    #   Diese Funktion stellt den größtmöglichen, "eindeutigsten" Eintrag 
    #   in der Liste S an vorderste Stelle. Sollte kein Eintrag eindeutig 
    #   sein, so wird der Eintrag der am seltensten vorkommt an vorderste
    #   Stelle gebracht. Ist auch das nicht eindeutig, so kommt der größte
    #   Wert an die erste Stelle

    count=0                 #Länge eines Blocks gleicher Zahlen
    num=S[0]                #Vorheriger Wert
    minCount=infinity
    indMinCount=0           #an welcher Stelle beginnt der kürzeste Block
    for i in range(0,len(S)):
        if S[i]==num:
            count=count+1
        else:   #d.h. ein Block ist zu Ende und der nächste hat begonnen
            if count<minCount:
                minCount=count
                indMinCount=i-count;
            if(count==1):   #kleinste Gruppe wurde schon jetzt sicher gefunden
                break;
            num=S[i]
            count=1

    #für den letzten Eintrag alles nochmal machen
    if count<minCount:
        minCount=count
        indMinCount=i+1-count
    #ein Element des kürzesten Blocks nach vorne schieben
    tmp=S[0]
    S[0]=S[indMinCount]
    S[indMinCount]=tmp
    return S

def Belyi_AGS(SW,SB,SI,info=False):
    r"""
    This function returns a system of equations, constructed with the signature
    of a dessin denfant, to calculate Belyi-functions. 

    AUTHORS:

    - Michael Eskin

    - Matthias Heinlein

    - Christian Steck

    INPUT:

    - ``SW`` a list with given ramificationindices of white points, which will
    be mapped to 0

    - ``SB`` a list with given ramificationindices of black points, which will
    be mapped to 1

    - ``SI`` a list with given ramificationindices of points, which will be
    mapped to infty

    - Optional: Bool "info"(=False by default). If info=True, all additional
    informations regarding normalization are        printed.

    OUTPUT:

    - A tuple ``[a,f]``, where ``a`` is a system of equations for the
    coefficients of the formal Belyi-function ``f``.

    NOTES: 

    - Our normalization is such, that f(0)=0, f(1)=1 and f(infty)=infty. Here
    the point x=0 is a white point with the 'most  unique', largest
    ramificationindex, and similar for x=1 and x=infty. If there is no unique
    value, then the highest value with the lowest crebritude will be mapped to
    the corresponding point.

    - To calculate ``f`` one has to solve the system of equations ``a`` and map
    these solutions to the coefficients of ``f``. (compare the function
            belyi-groebner). 

    .. WARNING::

        The lists SW,SB and SI should be entered with care. It is often
        advisable to use the list with the least elements or with the highest
        value as SI, since the degree of the denominator will be minimal. But
        in general, an automated change has proven to lead to bigger systems of
        equations. I.e. is it possible, that interchanging the lists can lead
        to a reduced need of computational power. The order of the elements in
        each list, however, does not matter.

    EXAMPLES:
     
    In this example we construct a Belyi-function with signature
    SW=[2,2],SB=[3,1],SI=[3,1]::

        sage: [a,f]=Belyi_AGS([2,2],[3,1],[3,1])
        sage: a
        [2*a0*a1 - a0*a3 + 3*a0, a0*a1^2 + 3*a0*a3 - 3*a0, a0*a3 - a2, -3*a0*a3 + a0 - 1]
        sage: b=solve_algebraic(a)
        [a0 + 3/8*a3 + 19/8, a1 - 1/2*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 + 3/8*a3 + 19/8, a1 + (-1/2)*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?, a2 + 0.3169872981077807?]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?]
        [a0 - 0.04903810567665797?]
        []
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?, a2 + 1.183012701892220?]
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?]
        [a0 + 2.549038105676658?]
        []

    This is the formal Belyi-function ``f``::

        sage: f
        (a0*x^4 + 2*a0*a1*x^3 + a0*a1^2*x^2)/(x + a2)

    We now map the coefficients of ``f`` to the solutions::

        sage: f=f.numerator().map_coefficients(b[0])/f.denominator().map_coefficients(b[0])
        sage: f
        (0.04903810567665797?*x^4 - 0.4641016151377546?*x^3 + 1.098076211353316?*x^2)/(x - 0.3169872981077807?)
    
    We test for the correctness of the ramification::

        sage: f.numerator().roots()
        [(0, 2), (4.732050807568877?, 2)]
        sage: (f-1).numerator().roots()
        [(1, 3), (6.464101615137755?, 1)]


    If "info=True" the first's command output differs, everything else remains
    the same::
        
        sage: [a,f]=Belyi_AGS([2,2],[3,1],[3,1],info=True)
        
        The chosen normalization:
        
        
        A white point with the rarest ramificationindex  2  is mapped to 0
        
        The unique black point with ramificationindex  3  is mapped to 1
        
        The unique point with ramificationindex  3  is mapped to infty 
        
        The resulting denominator of the Belyi-function has degree  1 
        
        The ring of the coefficients is 
        Multivariate Polynomial Ring in a0, a1, a2, a3 over Rational Field 
        
        The ring, where the calculations are actually done is 
        Univariate Polynomial Ring in x over Multivariate Polynomial Ring in a0, a1, a2, a3 over Rational Field
        
 
   """

    #Sortierung aller Listen um sie an eindToFront zu übergeben
    SW.sort(reverse=True)
    SB.sort(reverse=True)
    SI.sort(reverse=True)
    SW=eindToFront(SW)
    SI=eindToFront(SI)
    SB=eindToFront(SB)
    n=sum(SW)
    if( sum(SB)!=n or sum(SI) != n or sum(SB)!=sum(SI)):
        sys.exit(['Due to the fact, that the lists are no partition of ',n,' there cannot be a Belyi-morphism of this type. Please check your lists.']);
    #Riemann-Hurwitz
    if(-2!=(n-len(SW)-len(SB)-len(SI))):
        sys.exit('The lists contradict the formula of Riemann-Hurwitz.')


    if info==True:
        print '\nThe chosen normalization:\n\n'
        if len([o for o in SW if o==SW[0]])==1:
            print 'The unique white point with ramificationindex ',SW[0], ' is mapped to 0\n'
        else:
            print 'A white point with the rarest ramificationindex ',SW[0],' is mapped to 0\n'
        if len([o for o in SB if o==SB[0]])==1:
            print 'The unique black point with ramificationindex ',SB[0], ' is mapped to 1\n'
        else:
            print 'A black point with the rarest ramificationindex ',SB[0],' is mapped to 1\n'
        if len([o for o in SI if o==SI[0]])==1:
            print 'The unique point with ramificationindex ',SI[0], ' is mapped to infty \n'
        else:
            print 'A point with the rarest ramificationindex in SI, namely ',SI[0], 'is mapped to infty\n'
        print 'The resulting denominator of the Belyi-function has degree ', sum(SI)-SI[0],'\n'
   
    #Konstruktion der Polynomringe mit der richtigen Anzahl benötigter Variablen
    anzahl=len(SW)+len(SB)+len(SI)
    T=PolynomialRing(QQ,'a',anzahl-2,order='lex') #Order = 'lex' für Gröbner-Basis, 
                                                  #sonst, sehr funkioniert nicht.

    R.<x>=T[]

    if info==True:
        print 'The ring of the coefficients is \n',T,'\n'
        print 'The ring, where the calculations are actually done is \n',R,'\n'

    #Die gesuchte Funktion f hat die Gestalt f=k*A/C mit einer Konstanten k, die 
    #in der Realisierung von uns T.gen(0)=a0 ist.

    #Wegen der Sortierung haben wir einen optimalen Normalisierungsparameter
    A=x^SW[0]
    sw0=SW[0]   #Speichern des Wertes der gelöscht wird
    SW.remove(SW[0])


    #Benötigen zum erzeugen der Polynome A, C und nacher noch B die echt verschiedenen 
    #Verzweigungsindizes
    SWset=list(set(SW))
    SWset.sort(reverse=True)

    SBset=list(set(SB))
    SBset.sort(reverse=True)

    #Wir haben SI absteigend sortiert. Der höchste Verzweigungsindex wird auf infty geschickt.
    #Damit ist der Grad im Nenner minimal.
    SI.remove(SI[0])
    SIset=list(set(SI))
    SIset.sort(reverse=True)

    offset=0    #Welche Variablen benötigt man

    #Hier erstellen wir das Polynom A. Wir fassen gleiche Verzweigungsindizes zu Polynomen zusammen. 
    #z.B. [3,2,2] entspricht X^3*(X^2+aX+b)^2, [3,2,2,2] entspricht X^3(X^3+aX^2+bX+c)^2
    for a in SWset:
        d=len([o for o in SW if o==a])
        P=R(0)
        for i in range(d):
            P=P+T.gen(1+i+offset)*x^i
        offset=offset+d
        P=P+x^d
        P=P^a
        A=A*P
     
    #Analog für das Nennerpolynom C.
    C=R.one()
    if(len(SI)>=1):
        for a in SIset:
            d=len([o for o in SI if o==a])
            P=R(0)
            for i in range(d):
                P=P+T.gen(1+i+offset)*x^i
            offset=offset+d
            P=P+x^d
            P=P^a
            C=C*P

    #Um die Verzweigungen bei 1 festzulegen, faktorisieren wir die Nullstellen von f-1.
    #Aufgrund der Sortierung, nehmen wir den optimalen Verzweigungsindex und schicken den auf 1.
    B=(x-1)^SB[0]
    SB.remove(SB[0])

    #Hier erzeugen wir analog zu vorher das Polynom B.
    for a in SBset:
        d=len([o for o in SB if o==a])
        P=R(0)
        for i in range(d):
            P=P+T.gen(1+i+offset)*x^i
        offset=offset+d 
        P=P+x^d
        P=P^a
        B=B*P

    #Hier findet der Koeffizientenvergleich statt. f-1=k*A/C-1 == k*B/C.
    #Umgestellt ist dies also A-B == C. 
    a=koeffvergleich(T.gen(0)*A-T.gen(0)*B,C)
    #Es wird das Gleichungssystem 'a' und die symbolische Belyi-Funktion f=k*A/C ausgegeben.
    #Man muss nur noch die Lösungen des Gleichungssystems, also die Koeffizienten in f einsetzen
    #und man erhält die Belyi-Funktionen.
    return [a,T.gen(0)*A/C]

def koeffvergleich(f,g):
    #Diese Funktion führt einen Koeffizientenvergleich durch
    #und gibt das resultierende Gleichungssystem zurück
    if f.degree()==g.degree():
        return [f.coeffs()[i]-g.coeffs()[i] for i in range(len(f.coeffs())) if(f.coeffs()[i]-g.coeffs()[i]!=0)]
    else: #d.h. der Grad ist nicht gleich, dann fülle sinnvoll mit 0-en auf.
        if f.degree()<g.degree():
            h=g
            g=f
            f=h
        a=[]
        for i in range(1,f.degree()-g.degree()+1):
            a.append(f.coeffs()[-i])
        f=f.truncate(g.degree()+1)
        for i in range(len(f.coeffs())):
            if(f.coeffs()[i]-g.coeffs()[i]!=0):
                a.append(f.coeffs()[i]-g.coeffs()[i])
        return a


def belyi_groebner(SW,SB,SI,info=False):
    r"""
    This function uses Belyi_AGS(), to calculate Belyi-functions out of a given
    signature of a dessin denfant. 

    AUTHORS:

    - Michael Eskin

    - Matthias Heinlein

    - Christian Steck

    INPUT:

    - ``SW`` a list with given ramificationindices of white points, which will
    be mapped to 0

    - ``SB`` a list with given ramificationindices of black points, which will
    be mapped to 1

    - ``SI`` a list with given ramificationindices of points, which will be
    mapped to infty

    - Optional: Bool "info"(=False by default). If info=True, all additional
    informations regarding normalization are printed.

    OUTPUT:

    - A set ``M`` consisting out of Belyi-functions.

    NOTES: 

    - Our normalization is such, that f(0)=0, f(1)=1 and f(infty)=infty. Here
    the point x=0 is a white point with the 'most unique', largest
    ramificationindex, and similar for x=1 and x=infty. If there is no unique
    value, then the highest value with the lowest crebritude will be mapped to
    the corresponding point.

    EXAMPLES:

    This is a calculation of the Belyi-morphisms with signature
    SW=[2,2],SB=[3,1] and SI=[3,1]::

        sage: M=belyi_groebner([2,2],[3,1],[3,1])
        [a0 + 3/8*a3 + 19/8, a1 - 1/2*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 + 3/8*a3 + 19/8, a1 + (-1/2)*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?, a2 + 0.3169872981077807?]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?]
        [a0 - 0.04903810567665797?]
        []
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?, a2 + 1.183012701892220?]
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?]
        [a0 + 2.549038105676658?]
        []
        sage: M
        [(0.04903810567665797?*x^4 - 0.4641016151377546?*x^3 + 1.098076211353316?*x^2)/(x - 0.3169872981077807?), (-2.549038105676658?*x^4 + 6.464101615137755?*x^3 - 4.098076211353316?*x^2)/(x - 1.183012701892220?)]

    The same example with info=True::

        sage: M=belyi_groebner([2,2],[3,1],[3,1],info=True)

        The chosen normalization:


        A white point with the rarest ramificationindex  2  is mapped to 0

        The unique black point with ramificationindex  3  is mapped to 1

        The unique point with ramificationindex  3  is mapped to infty 
        
        The resulting denominator of the Belyi-function has degree  1 
        
        The ring of the coefficients is 
        Multivariate Polynomial Ring in a0, a1, a2, a3 over Rational Field 
        
        The ring, where the calculations are actually done is 
        Univariate Polynomial Ring in x over Multivariate Polynomial Ring in a0, a1, a2, a3 over Rational Field 
        
        [a0 + 3/8*a3 + 19/8, a1 - 1/2*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 + 3/8*a3 + 19/8, a1 + (-1/2)*a3 + 3/2, a2 + 1/8*a3 + 9/8, a3^2 + 6*a3 - 3]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?, a2 + 0.3169872981077807?]
        [a0 - 0.04903810567665797?, a1 + 4.732050807568877?]
        [a0 - 0.04903810567665797?]
        []
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?, a2 + 1.183012701892220?]
        [a0 + 2.549038105676658?, a1 + 1.267949192431123?]
        [a0 + 2.549038105676658?]
        []
        sage: M
        [(0.04903810567665797?*x^4 - 0.4641016151377546?*x^3 + 1.098076211353316?*x^2)/(x - 0.3169872981077807?), (-2.549038105676658?*x^4 + 6.464101615137755?*x^3 - 4.098076211353316?*x^2)/(x - 1.183012701892220?)]
        
    """
    #Das ist ein einfaches Programm, dass Belyi-Morphismen mit Hilfe des Mini-Projekts löst. 
    #Achtung, das ist NICHT sehr effektiv.
    [a,f]=Belyi_AGS(SW,SB,SI,info)
    b=solve_algebraic(a)
    morphisms=[]
    for i in range(len(b)):
        morphisms.append(f.numerator().map_coefficients(b[i])/f.denominator().map_coefficients(b[i]))
    return morphisms


def evalPower(F,f,e):
    """
    TESTS::

        sage: S.<x> = QQ[]
        sage: R.<y> = S[]
        sage: F = (x + 1)*y
        sage: f = x^2 + x + 1
        sage: evalPower(F,f,0)
        1
        sage: evalPower(F,f,1)
        2
        sage: evalPower(F,f,2)
        2
        sage: evalPower(F,f,3)
        1
        sage: F = (x^2+1)*y^3+y
        sage: evalPower(F,f,3)
        10
    """
    ret = 0
    power = f.parent().one()
    for v in F.coeffs():
        for (w,e1) in zip(v.coeffs(),range(v.degree()+1)):
            e2 = e-e1
            if (e2<0): continue
            ret += power[e2]*w
        power*=f
        power=power.truncate(e+1)
    return ret

#@func_persist
def newtonSchritt(F,fn,n,m):
    if m<=0: return fn
    R=F.parent() #C[x][y]
    S=R.base() #C[x]
    y=R.gen() 
    fnm=S(fn) #soll in C[x] liegen (auflösen von F nach y in Abhängigkeit von x!)
    abl=F.derivative(y)
    abl=abl(fnm(0),0) #Wir brauchen nur eine erste Approximation der Ableitung
    for i in range(m):
        fnm-=(evalPower(F,fnm,n+i+1)*S.gen()^(n+i+1))/abl
    return fnm

#@func_persist
def newtonSchritt2(F,fn,n,m):
    if m<=0: return fn
    R=F.parent() #C[y][x]
    S=R.base() #C[x]
    y=R.gen() 
    fnm=S(fn) #soll in C[x] liegen (auflösen von F nach y in Abhaengigkeit von x!)
    abl=F.derivative(y)
    abl=abl(fnm(0),0) #Wir brauchen nur eine erste Approximation der Ableitung
    for i in range(m):
        help=S(F(fnm)) #Damit das truncate richtig funktioniert muss man den Ring wechseln
        fnm=fnm-help.truncate(n+i+2)/abl #Das truncate verhindert, dass die Lösung explodiert!
    return fnm
    
def newtonSchritt3(F,fn,n,m):
    if m<=0: return fn
    R=F.parent() #C[y][x]
    S=R.base() #C[x]
    x=S.gen()
    y=R.gen() 
    fnm=S(fn) #soll in C[x] liegen (auflösen von F nach y in Abhaengigkeit von x!)
    abl=F.derivative(y)
    T=PowerSeriesRing(S.base(),S.gens())
    i=n+1
    while i<m+n+1:
        print(i)
        i=2*i
        help=S(F(fnm)) #Damit das truncate richtig funktioniert muss man den Ring wechseln
        abl=T(abl(fnm)) #Wir brauchen nur eine erste Approximation der Ableitung
        abl=(~abl).truncate(i)
        fnm=fnm-help*abl #Das truncate verhindert, dass die Lösung explodiert!
    return fnm

def wurzelPolynom(P,e,m,k=-1,Base=QQbar):
    if e<=1:
        return P     
    
    S=Base['x']
    x=S.gen()
    P=S(P)
    if(k<0):k=P.valuation()
    P=P.shift(-k)            #Von x kann man keine Wurzel ziehen. Diese Terme werden spaeter dazugefügt
    R=S['y']
    y=R.gen()
    F=y^e-P
    
    P0=Base(P(0)^(1/e))
    
    #if k>0: P0=x^(k/e)*P0
    
    return newtonSchritt(F,P0,0,m-k/e)*x^(k/e)
    
    
def umkehrung(A,B,x0,y0,e,m,n,Base=CC):
    R=QQbar['x']
    A=R(A)
    B=R(B)
    x=R.gen()
    B=R(B(x+x0))
    A=R(A(x+x0))#nun ist das Problem auf x0=0 zurückgeführt
    A=R(A-y0*B)#nun ist das Problem auf x0=0 und y0=0 zurückgeführt
    
    A=wurzelPolynom(A,e,m,k=e,Base=Base)
    B=wurzelPolynom(B,e,m,k=0,Base=Base)
    
    
    R=A.parent()
    S=Base['y']['x']
    xneu=S.gen()
    y=S.base().gen() 
    phi=R.hom([xneu])
    #nach der richtigen Variable auflösen!

   
    P=phi(A)-y*phi(B) #Gleichung, die wir lösen wollen
            #A(0)=0 aber B(0)!=0
            
    P=newtonSchritt(P,0,0,n)
    return P+x0 #Das ist jetzt in Abhaengigkeit vom Parameter y_new, welcher etwa (y-y0)^(1/e) ist! Es gilt also (approx) A(P)/B(P)=y^e+A(x0)/B(x0).
    
def wertAn(P,y0,y,e):
    y=CC((y-y0)^(1/e))
    l=Einheitswurzeln(e)   
    return [P(h*y) for h in l]
    
def Einheitswurzeln(e):
    #S.<x>=QQbar[]
    #f=x^e-1
    #return [h[0] for h in f.roots()]    
    return[CC(exp(2*pi*i*k/e)) for k in range(e)]

def Zuordnung(P,r,y0,e,Pliste):
    if((y0<1/2)&(y0+r>1/2)):
        r=1/2-y0
    if((y0>1/2)&(y0+r<1/2)):
        r=1/2-y0    
    l=wertAn(P,y0,y0+r,e)
    l2=[Q(y0+r-1/2) for Q in Pliste]
    ret=[]
    for wert in l:
        opt=abs(wert-l2[0])
        iopt=0
        for i in range(len(l2)):
            if abs(wert-l2[i])<=opt:
                opt=abs(wert-l2[i])
                iopt=i
        ret.append(iopt)
    return ret
    
def konvergenzradius(A,B,P,y0,r,e,eps=0.001,abr=1/100):
    fehler=eps+1
    r=r*10/9
    while(fehler>eps):
        r=r/10*9
        fehler=max([abs(CC(A(h)/CC(B(h)))-y0-r) for h in wertAn(P,y0,y0+r,e)])
        #if(r<abr):
        #    return -1
    end
    return CC(r)    
    
    
def graph(A,B,eps=0.001,n0=10,n1=10,m=50,abr=1/100,Base=CC,komment=False):
    R.<x>=QQbar[]
    A=R(A)
    B=R(B)
    C=R.gcd(A,B)
    A=R(A/C)
    B=R(B/C)
    H=A-B
    l=H.roots()
    V1=[h[0] for h in l]
    e1=[h[1] for h in l]
    l=A.roots()
    V0=[h[0] for h in l]
    e0=[h[1] for h in l]
    H=A-1/2*B
    l=H.roots()
    l=[h[0] for h in l]
    e12=[1 for h in l]
    
    if(komment):print('Umkehrfunktion an 1/2 wird berechnet')
    Pliste=[umkehrung(A,B,l[i],1/2,e12[i],m,m,Base=Base) for i in range(len(l))] #Umkehrfunktionen für 1/2
    if(komment):print('Umkehrfunktion an 1 wird berechnet')
    P1=[umkehrung(A,B,V1[i],1,e1[i],n1,n1,Base=Base) for i in range(len(V1))] #Umkehrfunktionen für 1
    if(komment):print('Umkehrfunktion an 0 wird berechnet')
    P0=[umkehrung(A,B,V0[i],0,e0[i],n0,n0,Base=Base) for i in range(len(V0))] #Umkehrfunktionen für 0
    
    
    if(komment):print('Konvergenzradien schaetzen')    
    r1=[konvergenzradius(A,B,P1[i],1,-1,e1[i],eps=eps,abr=abr) for i in range(len(P1))]
    r0=[konvergenzradius(A,B,P0[i],0,1,e0[i],eps=eps,abr=abr) for i in range(len(P0))]
    r2=[konvergenzradius(A,B,h,1/2,1,1,eps=eps,abr=abr) for h in Pliste]
    if(komment):
        print('kleinster Konvergenzradius bei 1/2:')
        print(min(r2))
        print('kleinster Konvergenzradius bei 1:')
        print(-max(r1))
        print('kleinster Konvergenzradius bei 0:')
        print(min(r0))
        
        print('Werte zuordnen')
    
    l1=[]
    l0=[]
    for i in range(len(P1)):
        l1.extend(Zuordnung(P1[i],r1[i],1,e1[i],Pliste))
    for i in range(len(P0)):
        l0.extend(Zuordnung(P0[i],r0[i],0,e0[i],Pliste)) 
    V1neu=[]
    V0neu=[]
    for i in range(len(V1)) :
        for j in range(e1[i]):
            V1neu.append(i)
    for i in range(len(V0)):
        for j in range(e0[i]):
            V0neu.append(i)
    
    sigma0=[]
    sigma1=[]
    for i in range(len(V0)):
        sigma0.append(zyklus(V0neu,V1neu,l0,l1,i,0))
    for i in range(len(V1)):
        sigma1.append(zyklus(V0neu,V1neu,l0,l1,i,1))
    
      
     
    if(min(r2)-max(r1)<1/2 or min(r2)+min(r0)<1/2):
        print('DIE GENAUIGKEIT REICHT NICHT AUS, UM EINE ZUORDNUNG MIT DEN ANGEGEBEN FEHLERN ZU MACHEN!!!!!')
        print('Man veraendere dazu den optionalen Parameter n0, n1, m (Genauigkeit fuer die Konvergenzreihen um 0 und 1 und 1/2, default 10,10,50).')
        print('Bei kleinen Konvergenzradien kann dies auch an der verwendeten Genauigkeit liegen. Benutze optionales Argument Base=CDF.to_prec(...)!')
        print('Minimaler Konvergenzradius:')
        print(min([min(r2),-max(r1),min(r0)]))

     
    return [sigma0,sigma1,V0,V1,l,A,B]
    
def zyklus(V0,V1,l0,l1,i,k):
    if k!=1: 
        return [[h[1],h[0],h[2]] for h in zyklus(V1,V0,l1,l0,i,1)]
    ret=[]
    for l in range(len(V1)):
        if V1[l]==i:
            for w in range(len(l0)):
                if l1[l]==l0[w]:
                    kante=[V0[w],i,l1[l]]
                    ret.append(kante)
    return ret
    
def Kinderzeichnung(liste,n=1,m=50,Base=CC,eps=0.1,r=-1,rmax=0.1):
    g = Graphics()
    A=liste[5]
    B=liste[6]
    V0=liste[2]
    V=liste[3]
    V.extend(V0)
    if(r<0):
        r=min([rmax,min([help for help in [abs(CC(V[i]-V[j])) for i in range(len(V)) for j in range(len(V))] if (not help==0)])/3])
    for k in getkanten(liste):
        g=gutekante(g,A,B,liste[2][k[0]],liste[4][k[2]],liste[3][k[1]],n,m,r,Base=Base,eps=eps)    
    g.axes(False)
    return g

def kante(g,path,r):
    path=[CC(p) for p in path]
    path=[(p.real_part(),p.imag_part()) for p in path]
    g+=line(path, linestyle='solid', rgbcolor='black')
    g += circle(path[0],r, fill=True, facecolor='white',edgecolor='black')
    g += circle(path[len(path)-1],r,color='black', fill=True)          
    return g
    
def gutekante(g,A,B,v0,l,v1,n,m,r,Base=CC,eps=0.001):
    path=[v0]
    if(n>1):
        y=[0,1/2,1]
        while(len(y)<n+2):
            s=len(y)
            for i in range(s):
                if(not i==s-1):
                    j=s-i-1
                    y.insert(j,(y[j]+y[j-1])/2)
        P=umkehrung(A,B,l,1/2,1,m,m,Base=Base)
        y.pop()
        y.remove(0)
        x=[wertAn(P,1/2,y0,1)[0] for y0 in y]
        for i in range(len(x)):    
            if(abs(A(x[i])/B(x[i])-y[i])>eps):
                print("In der Zeichnung wurden Punkte benutzt, die nicht gut approximiert wurden!!!")
                print("optinale Parameter dafuer: n (Zahl der Zwischenpunkte - klein gut, dafuer aber eventuell Kreuzungen) m (Approximationsguete der Umkehrfunktion)")
        path.extend(x)
        path.append(v1)
    if(n==1):
        path=[v0,l,v1]
    g=kante(g,path,r)
    return g
    
def getkanten(liste):
    ret=[]
    for z in liste[0]:
        for k in z:
            if(k not in ret): 
                ret.append(k)
    for z in liste[1] :
        for k in z :
            if(k not in ret): 
                ret.append(k)
    return ret
    
def getsigma(i,liste):
    help=[[k[2] for k in h] for h in liste[i]]
    ret=Permutation([])
    for h in help:
        ret=ret*getPermutation(h)
    return ret
        
def getsigmaInf(liste):
    s1=getsigma(1,liste)
    s0=getsigma(0,liste)
    s=s0*s1
    return s.inverse()

def getPermutation(h):
    ret=[f+1 for f in range(max(h)+1)]
    for i in range(len(h)):
        if(i<(len(h)-1)):ret[h[i]]=h[i+1]+1
        if(i==(len(h)-1)):ret[h[i]]=h[0]+1
    return(Permutation(ret))


#Hier hängt noch das Miniprojekt dran.
def solve_algebraic(equations):
    R=equations[0].parent()
#    assert(equations[0].base_ring()==QQ)
    I=R.ideal(equations)
    G=I.groebner_basis()
    print G
#Lift auf QQbar
    RQbar=PolynomialRing(QQbar,list(R.gens()))
    liftqqbar=hom(R,RQbar,list(RQbar.gens()))
    Gneu=[1]*len(G)
    for i in range(0,len(G)):
        Gneu[i]=liftqqbar(G[i])
    return solve_algebraic_groebner(Gneu)

def solve_algebraic_groebner(G):
    print G
    if len(G) == 0:
        return [QQbar.hom(QQbar)]
    else:
        ret = []
        g=G[-1]
        S= solve_algebraic_univariate(make_univariate(g))
        for s in S:
            s=lift_with_identity(s,g.parent())
            blub=solve_algebraic_groebner([s(h) for h in G[0:-1]])
            blub=[fun*s for fun in blub] 
            ret=ret+blub
        return ret

def make_univariate(p):
    """
    EXAMPLES::

        sage: R.<x> = QQbar[]
        sage: f = x
        sage: g = make_univariate(f)
        sage: g
        x
        sage: g.parent()
        Univariate Polynomial Ring in x over Algebraic Field

        sage: R.<x,y> = QQbar[]
        sage: f = y+1
        sage: g = make_univariate(f)
        sage: g
        y + 1
        sage: g.parent()
        Univariate Polynomial Ring in y over Algebraic Field
    """
    return p if p.parent().ngens()==1 else p.univariate_polynomial()

def lift_with_identity(s,domain):

#neuen Ring erzeugen
    anzahl=domain.ngens()-s.domain().ngens()
    neuevariablen=list(domain.gens())[0:anzahl]
    Rneu=PolynomialRing(QQbar,neuevariablen)

    insert=list(Rneu.gens())+s.im_gens()

    f=domain.hom(insert)

    return f

def solve_algebraic_univariate(equation):

# assert(equation.parent().ngens()==1)
#    assert(not equation.is_constant())

    return [equation.parent().hom([r]) for r in equation.roots(multiplicities=false)]
