def p_adic_newton(p,f,n):
    r"""
    Returns p-adic solutions of a set of polynomial equations f

    DESCRIPTION:
    
    This function first calculates solutions of a set of integer polynomial 
    equations f modulo p. Using Newton's Algorithm we obtain a p-adic approximative
    solution.

    AUTHORS:

    - Matthias Rapp: initial version

    - Michael Zell: 


    EXAMPLES:

2 Equations::
    
    sage: Q.<x,y>=QQ[]
    sage: f=range(2)
    sage: f[0]=x^2-2
    sage: f[1]=y^2-4
    sage: p_adic_newton(7,f,10)
    
::
3 Equations::

    sage: Q.<x,y,z>=QQ[]
    sage: f=range(3)
    sage: f[0]=x^2-y+1
    sage: f[1]=x^2-y+z
    sage: f[2]=1-y+z
    sage: p_adic_newton(7,f,10)

::
    """
    if not(isinstance(f,list)):
        f=[f]
        
    Q=f[0].parent()
    if len(f)!=len(Q.gens()):
        print "ERROR: Number of equations and variables doesn't match"
        return   
   
    
    Z=Zp(p) 
    g=range(len(f))
    for i in range(len(f)):
        g[i]=f[i].change_ring(Z)
    
    R=Q.change_ring(Z)
    J=matrix(R,len(f))
    K=matrix(R,len(f),1)
    
    # Create Matrix for substitution
    for i in range(len(f)):
        K[(i,0)]=g[i]    
    
    # Calculate Jacobian Matrix
    for i in range(len(f)):
        for j in range(len(Q.gens())):
            h=diff(f[i],Q.gens()[j])   
            J[(i,j)]=h.change_ring(Z) 

    # Calculate solutions mod p
    A=Q.change_ring(GF(p))
    I=A.ideal(f)
    V=I.variety()
    sol=[[v[t] for t in A.gens()] for v in V]
    
    lst=[] # Storage for Results
    
    # Perform Newton Algorithm
    for h in range(len(sol)):
        a=[]
        for i in range(len(Q.gens())):
            a.append(Z(sol[h][i].lift()))
        
        for k in range(n-1):
            if GF(p)(det(J)(a)) != 0:
                delta = ~J(a)*(-K(a))
                tmp=matrix(len(f),1,a)+delta
                a=tmp.list()
            else:
                print "No Convergence possible, Jacobian Matrix is singular"
                return
        lst.append(a)
    # Return list of Results
    return lst



# Calculate minimal Polynom for x using LLL-Algorithm
def mipo(x,n):
    r"""
    Returns the minimal polynoms for a set x of p-adic numbers

    DESCRIPTION:
    
    Calculates the minimal polynoms for a set x of p-adic numbers
    using the LLL Algorithm. 
    
    n: maximal degree the minimal polynom can reach
    power: large number, e.g. if you use Zp() => power=20 is the accuracy of Zp()
 

    AUTHORS:

    - Matthias Rapp: initial version

    - Michael Zell: 


    EXAMPLES:

Minimal polynom for p-adic construction of 3/4::

    sage: Z=Zp(7)
    sage: mipo(Z(3/4),4)
    [4*X - 3]  

::

    """
    if not(isinstance(x,list)):
        x=[x]
        p=x.parent().base_ring().prime()
        power=x.parent().base_ring().precision_cap()
    else:
        elt=x[0]
        p=elt.parent().base_ring().prime()
        power=elt.parent().base_ring().precision_cap()

    M=Matrix(ZZ,n,n+1)
    Q.<X>=ZZ[]
    ans=[]
    
    # Perform LLL on every component
    for j in range(len(x)):
        for i in range(n):
            M[i,i]=1

        for i in range(n):
            M[i,n]=p^power*x[j]^(i+1)

        L=M.LLL()

        P=0
    
        for i in range(n):
            P+=X^i*L[0,i]
        
        # Factor candidate for minimal polynom
        fac=list(P.factor())
        # Look for factor which has appropriate root
        for i in range(len(fac)):
            k=fac[i][0](x[j])
            if ZZ(k.lift())==0:
                 ans.append(fac[i][0])
    return ans
    
def solve_p_adic(p,f,n):
    r"""
    Returns homomorphisms as solutions for an set of equations f

    DESCRIPTION:
    
    See above.
    

    AUTHORS:

    - Matthias Rapp: initial version

    - Michael Zell: 


    EXAMPLES:

Minimal polynom for p-adic construction of 3/4::

    sage: Q.<x,y>=QQ[]
    sage: f=range(2)
    sage: f[0]=x^2-2
    sage: f[1]=y^2-4
    
    sage: solve_p_adic(7,f,10)
        [Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Algebraic Field
          To:   Algebraic Field
          Defn: x |--> -1.414213562373095?
                y |--> 2, Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Algebraic Field
          To:   Algebraic Field
          Defn: x |--> 1.414213562373095?
                y |--> 2, Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Algebraic Field
          To:   Algebraic Field
          Defn: x |--> -1.414213562373095?
                y |--> -2, Ring morphism:
          From: Multivariate Polynomial Ring in x, y over Algebraic Field
          To:   Algebraic Field
          Defn: x |--> 1.414213562373095?
                y |--> -2]
::

    """
    solutions=p_adic_newton(p,f,n)
    tmp=[]  
    rootpairs=[]
    bla=[]
    for i in range(len(solutions)):
        tmp=[]
        mipos=mipo(solutions[i],n)
        for j in range(len(mipos)):
            tmp.append(mipos[j].change_ring(QQbar).roots())  
            
            a=range(len(tmp))
            for l in range(len(tmp)):
                a[l]=[]
                for k in range(len(tmp[l])):
                    a[l].append(tmp[l][k][0])
            bla=a
           
            del(a)
        for l in range(len(list(cartesian_product_iterator(bla)))):
            rootpairs.append(list(cartesian_product_iterator(bla))[l])  
    homs=[]
    for i in range(len(rootpairs)):
        homs.append(f[0].parent().change_ring(QQbar).hom(rootpairs[i],QQbar))
    return [x for i,x in enumerate(homs) if x not in homs[i+1:]]
