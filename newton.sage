#Q.<x1,x2>=QQ[]
#f=range(2)
#f[0]=x1^2-2
#f[1]=x2^2-4

def newton(Q,p,f,n):
    
    Z=Zp(p) 
    g=range(len(f))
    for i in range(len(f)):
        g[i]=f[i].change_ring(Z)
    
    R.<x1,x2>=Z[]
    J=matrix(R,len(f))
    K=matrix(R,len(f),1)

    # Create Matrix for substitution
    for i in range(len(f)):
        K[(i,0)]=g[i]    
    
    # Calculate Jacobian Matrix
    for i in range(len(f)):
        for j in range(len(f)):
            h=diff(f[i],Q.gens()[j])   
            J[(i,j)]=h.change_ring(Z) 
            
    # Calculate Solution mod 7     
    x=var(Q.variable_names())
    b=[]
    for i in range(len(f)):
        b.append(f[i](x))
    sol=solve_mod(b,p)
    
    lst=[] # Storage for Results
    
    # Perform Newton Algorithm
    for h in range(len(sol)):
        a=[]
        for i in range(len(Q.gens())):
            a.append(QQ(sol[h][i]))
    
        for k in range(n-1):
            if GF(p)(det(J(a))) != 0:
                delta = ~J(a)*(-K(a))
                tmp=matrix(len(f),1,a)+delta
                a=tmp.list()
            else:
                print "ERROR"
                break
        lst.append(a)
    # Return list of Results
    return lst

# Calculate minimal Polynom for x using LLL-Algorithm
def mipo(x,p,n,power):
    M=Matrix(ZZ,n,n+1)
    ans=[]
    for j in range(len(x)):
        for i in range(n):
            M[i,i]=1

        for i in range(n):
            M[i,n]=p^power*x[j]^(i)

        L=M.LLL()

        P=0
        X=var("X")
    
        for i in range(n):
            P+=X^i*L[0,i]
        ans.append(P)
    return ans

