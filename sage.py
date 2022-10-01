#!/usr/bin/env python
# coding: utf-8

# In[1]:


# from sage.all import *
import time
from sage.plot.plot3d.shapes import *

get_ipython().run_line_magic('display', 'latex')


# Definimos las coordenadas $(t=x^0,x^1,x^2,x^3)$ en que estará definida la métrica $g_{\mu\nu}$.
# 
# Definimos también la variable $v$ que codificará los órdenes de magnitud en velocidad de las componentes de la métrica y objetos que deriven de ella.
# 
# La función truncate coge un polinomio en $x$ y lo corta a orden $n$.
# 
# La función diff\_v deriva una funcion $f(t,x^1,x^2,x^3)$ respecto de una coordenada. Cuando se deriva respecto de $x^0=t$, dado que $\partial/\partial t\sim v/r$, se añade una unidad al orden de magnitud en v.

# In[2]:


inicio = time.time()

show = False    #True para display de la métrica, conexión afin, etc

metrica = True   #True para introducir una expansión en serie de una métrica y realizar los cálculos con ella

t_name='t'
x1_name='r'
x2_name='\\theta'
x3_name='\\phi'

v = var('v')
t = var('t', latex_name=t_name)
x1 = var('x1', latex_name=x1_name)
x2 = var('x2', latex_name=x2_name)
x3 = var('x3', latex_name=x3_name)

G, M, c = var('G M c')

n = 4 # Orden de magnitud en v al que queremos llegar

if metrica==True: #Introducir una métrica como un polinomio en v
    for i in range(4):
        for j in range(i,4):
            exec('in_g'+str(i)+str(j)+' = function("in_g"+str(i)+str(j))(t,x1,x2,x3)')
            
    in_g00 = -c**2+2*G*M/x1*v**2
    in_g01 = 0*v
    in_g02 = 0*v
    in_g03 = 0*v
    in_g11 = 1 + 2*G*M/c**2/x1*v**2 + (2*G*M/c**2/x1)**2*v**4
    in_g12 = 0*v
    in_g13 = 0*v
    in_g22 = x1**2
    in_g23 = 0*v
    in_g33 = x1**2*sin(x2)**2
    
    #variables para los terminos de orden 0 en la diagonal
    if in_g00.coefficients(v)[0][1]==0:
        orden0_0 = in_g00.coefficients(v)[0][0]
    if in_g11.coefficients(v)[0][1]==0:
        orden0_1 = in_g11.coefficients(v)[0][0]
    if in_g22.coefficients(v)[0][1]==0:
        orden0_2 = in_g22.coefficients(v)[0][0]
    if in_g33.coefficients(v)[0][1]==0:
        orden0_3 = in_g33.coefficients(v)[0][0]

else:
    orden0_0 = -1
    orden0_1 = 1
    orden0_2 = 1
    orden0_3 = 1

aux = var('aux') #Variable auxiliar para poder definir la matriz gseries como una matriz de funciones y no de números, luego
                    #se evalua en aux=0
    
def truncate(pol, x, n): #Función para truncar un polinomio pol(x) a orden n
    i=0
    pol_trunc=0
    if(len(pol.coefficients(x))==0):
        return pol
    if(pol.coefficients(x)[len(pol.coefficients(x))-1][1]>n):
        while(pol.coefficients(x)[i][1] <= n):
            i+=1
    else:
        i=len(pol.coefficients(x))
    for j in range(i):
        pol_trunc += pol.coefficients(x)[j][0]*x**(pol.coefficients(x)[j][1])
    return pol_trunc

def diff_v(f, k): #Función para derivar f respecto de la coordenada k-ésima
    variables = [t, x1, x2, x3]
    der = function('der')(t,x1,x2,x3)
    if(k==0):
        der = diff(f, t, 1)*v
    else:
        der = diff(f, variables[k], 1)
    return der


# Definimos las componentes de la métrica como una expansión en serie de potencias (en orden de magnitud) de la velocidad,
# $$g_{00}=-1+g_{00}^{(2)}+g_{00}^{(4)}+...$$
# $$g_{ij}=\delta_{ij}+g_{ij}^{(2)}+g_{ij}^{(4)}+...$$
# $$g_{0i}=g_{0i}^{(3)}+g_{0i}^{(5)}+...$$

# In[3]:


# Definimos los coeficientes del desarrollo en serie de la métrica

if(n%2==0):
    m=n/2
else:
    m=(n+1)/2

if(n%2==0):  # Cuando n es par, los elementos 00 e ij tienen que correr hasta n/2+1 y los 0i hasta n/2. Cuando n es impar, todos hasta (n+1)/2
    m+=1
    
for k in range(1, m):
    exec('g00_'+str(2*k)+' = var("g00_"+str(2*k), latex_name="g^{("+str(2*k)+")}_{00}")')
    exec('f_g00_'+str(2*k)+' = function("f_g00_"+str(2*k), latex_name="g^{("+str(2*k)+")}_{00}")(t,x1,x2,x3)')

for i in range(1,4):
    for j in range(i,4):
        for k in range(1, m):
            exec('g'+str(i)+str(j)+'_'+str(2*k)+' = var("g"+str(i)+str(j)+"_"+str(2*k), latex_name="g^{("+str(2*k)+")}_{"+str(i)+str(j)+"}")')
            exec('f_g'+str(i)+str(j)+'_'+str(2*k)+' = function("f_g"+str(i)+str(j)+"_"+str(2*k), latex_name="g^{("+str(2*k)+")}_{"+str(i)+str(j)+"}")(t,x1,x2,x3)')

if(n%2==0):
    m-=1

for i in range(1,4):
    for k in range(1, m):
        exec('g0'+str(i)+'_'+str(2*k+1)+' = var("g0_"+str(i)+str(2*k+1), latex_name="g^{("+str(2*k+1)+")}_{0"+str(i)+"}")')
        exec('f_g0'+str(i)+'_'+str(2*k+1)+' = function("f_g0_"+str(i)+str(2*k+1), latex_name="g^{("+str(2*k+1)+")}_{0"+str(i)+"}")(t,x1,x2,x3)')


# In[4]:


gseries = Matrix(4,aux)

if(n%2==0):
    m+=1

# Definimos g00 en serie
gseries[0,0] += orden0_0    ################## += -1
for k in range(1, m):
    exec('gseries[0,0] += g00_'+str(2*k)+'*v**(2*k)')

# Definimos gij en serie
for i in range(1,4):
    for j in range(i,4):
        if(i==j):
            exec('gseries[i,i] += orden0_'+str(i))     ################## += 1
        for k in range(1, m):
            exec('gseries[i,j] += g'+str(i)+str(j)+'_'+str(2*k)+'*v**(2*k)')
        gseries[j,i] = gseries[i,j]
        
if(n%2==0):
    m-=1

# Definimos g0i en serie
for i in range(1,4):
    for k in range(1, m):
        exec('gseries[i,0] += g0'+str(i)+'_'+str(2*k+1)+'*v**(2*k+1)')
    gseries[0,i] = gseries[i,0]
            
# Eliminamos aux evaluando en 0
gseries = gseries(aux=0)
if show==True:
    display(gseries)


# Análogamente, definimos los coeficientes de la inversa de la métrica de la misma forma:
# $$g^{00}=-1+g^{00}_{(2)}+g^{00}_{(4)}$$
# $$g^{ij}=\delta_{ij}+g^{ij}_{(2)}+g^{ij}_{(4)}$$
# $$g^{0i}=g^{0i}_{(3)}+g^{0i}_{(5)}$$

# In[5]:


# Definimos los coeficientes del desarrollo en serie de la inversa de la métrica

if(n%2==0):
    m+=1

for k in range(1, m):
    exec('ginv00_'+str(2*k)+' = var("ginv00_"+str(2*k), latex_name="g^{00}_{("+str(2*k)+")}")')

for i in range(1,4):
    for j in range(i,4):
        for k in range(1, m):
            exec('ginv'+str(i)+str(j)+'_'+str(2*k)+' = var("ginv"+str(i)+str(j)+"_"+str(2*k), latex_name="g^{"+str(i)+str(j)+"}_{("+str(2*k)+")}")')

if(n%2==0):
    m-=1
            
for i in range(1,4):
    for k in range(1, m):
        exec('ginv0'+str(i)+'_'+str(2*k+1)+' = var("ginv0"+str(i)+"_"+str(2*k+1), latex_name="g^{0"+str(i)+"}_{("+str(2*k+1)+")}")')


# In[6]:


ginvseries = Matrix(4,aux)

if(n%2==0):
    m+=1

# Definimos g^00 en serie
ginvseries[0,0] += 1/orden0_0     ################## += -1
for k in range(1, m):
    exec('ginvseries[0,0] += ginv00_'+str(2*k)+'*v**(2*k)')

# Definimos g^ij en serie
for i in range(1,4):
    for j in range(i,4):
        if(i==j):
            exec('ginvseries[i,i] += 1/orden0_'+str(i))     ################## += 1
        for k in range(1, m):
            exec('ginvseries[i,j] += ginv'+str(i)+str(j)+'_'+str(2*k)+'*v**(2*k)')
        ginvseries[j,i] = ginvseries[i,j]

if(n%2==0):
    m-=1
        
# Definimos g^i0 en serie
for i in range(1,4):
    for k in range(1, m):
        exec('ginvseries[i,0] += ginv0'+str(i)+'_'+str(2*k+1)+'*v**(2*k+1)')
    ginvseries[0,i] = ginvseries[i,0]
            
# Eliminamos aux evaluando en 0
ginvseries = ginvseries(aux=0)
if show==True:
    display(ginvseries)


# Para expresar los coeficientes de la inversa en función de los $g_{ij}^{(k)}$ resolvemos el sistema de ecuaciones
# $$g^{\mu\rho}g_{\rho\nu}=\delta_{\mu\nu}$$

# In[7]:


# Calculamos el producto de la métrica y su inversa, truncamos todos los elementos a n-ésimo orden y resolvemos las ecuaciones
    # igualando a la matriz identidad, resolviendo para los coeficientes de la inversa

f=ginvseries*gseries
ginv=Matrix(4,aux)

for i in range(4):
    for j in range(4):
        f[i,j] = truncate(f[i,j], v, n)
        
    
# Recorremos la matriz f resolviendo cada ecuación desde órdenes menores hacia órdenes mayores y sustituyendo las soluciones

for k in range(1,n+1):
    for i in range(4):
        for j in range(i,4):
            long = len(f[i,j].coefficients(v))
            l = 0
            while(l < long):
                if(f[i,j].coefficients(v)[l][1] == k):
                    exec("solution=solve(f[i,j].coefficients(v)[l][0]==0, ginv"+str(i)+str(j)+"_"+str(f[i,j].coefficients(v)[l][1])+")")
                    if(show==True):
                        display(solution)
                    for q in range(4):
                        for p in range(4):
                            f[q,p] = f[q,p].subs(solution)
                            ginvseries[q,p] = ginvseries[q,p].subs(solution)
                    longnew = len(f[i,j].coefficients(v))
                    if(long != longnew):
                        l = l-1
                        long = longnew
                l = l+1


# Dado que para resolver el sistema anterior hemos tenido que definir los coeficientes de la métrica y la inversa como variables, y dado que a continuación tendremos que calcular sus derivadas, los sustituimos por funciones de las coordenadas $(t,x^1,x^2,x^3)$:

# In[8]:


#Si introducimos una métrica, sustituimos todas las funciones g_ij(k) por su valor introducido
#Si no introducimos una métrica, queda en función del desarrollo en serie calculado

g = Matrix(4,aux)
ginv = Matrix(4,aux)

for i in range(4):
    for j in range(i,4):
        exec('g[i,j] = g[j,i] = gseries[i,j]')
        exec('ginv[i,j] = ginv[j,i] = ginvseries[i,j]')
        

if metrica==True: #f_gij_k
    for i in range(4):
        for j in range(i,4):
            exec('coef_metrica = in_g'+str(i)+str(j)+'.coefficients(v)')
            coef=[0 for p in range(n)]
            if (i==0 and j!=0):
                if coef_metrica[0][0]==0:
                    for k in range(m):
                        coef[k]=[0,2*k+1]
                else:
                    alpha=0
                    for k in range(m):
                        if alpha<len(coef_metrica) and coef_metrica[alpha][1]==2*k+1:
                            coef[k]=coef_metrica[alpha]
                            alpha+=1
                        else:
                            coef[k]=[0,2*k+1]
            else:
                if(n%2==0):
                    m+=1
                if coef_metrica[0][0]==0:
                    for k in range(m):
                        coef[k]=[0,2*k]
                else:
                    alpha=0
                    for k in range(m):
                        if alpha<len(coef_metrica) and coef_metrica[alpha][1]==2*k:
                            coef[k]=coef_metrica[alpha]
                            alpha+=1
                        else:
                            coef[k]=[0,2*k]
                if(n%2==0):
                    m-=1
                
            if (i==0 and j!=0):
                for k in range(1,m):
                    exec('f_g'+str(i)+str(j)+'_'+str(2*k+1)+' = coef[k][0] ')
            else:
                if(n%2==0):
                    m+=1
                for k in range(1,m):
                    exec('f_g'+str(i)+str(j)+'_'+str(2*k)+' = coef[k][0] ')
                if(n%2==0):
                    m-=1
                


# In[9]:


# Reemplazamos las variables gij_k por funciones gij_k(t,x1,x2,x3) para poder derivarlas cuando sea necesario
          
if(n%2==0):
    m+=1

# replace 00
for i in range(4):
    for j in range(4):
        for k in range(1,m):
            exec('g[i,j] = g[i,j].subs( g00_'+str(2*k)+' == f_g00_'+str(2*k)+')')
            exec('ginv[i,j] = ginv[i,j].subs( g00_'+str(2*k)+' == f_g00_'+str(2*k)+')')
                
# replace ij
for i in range(4):
    for j in range(4):
        for k in range(1,m):
            for q in range(1,4):
                for p in range(q,4):
                    exec('g[i,j] = g[i,j].subs( g'+str(q)+str(p)+'_'+str(2*k)+' == f_g'+str(q)+str(p)+'_'+str(2*k)+')')
                    exec('ginv[i,j] = ginv[i,j].subs( g'+str(q)+str(p)+'_'+str(2*k)+' == f_g'+str(q)+str(p)+'_'+str(2*k)+')')

if(n%2==0):
    m-=1
                    
# replace 0i
for i in range(4):
    for j in range(4):
        for k in range(1,m):
            for p in range(1,4):
                exec('g[i,j] = g[i,j].subs( g0'+str(p)+'_'+str(2*k+1)+' == f_g0'+str(p)+'_'+str(2*k+1)+')')
                exec('ginv[i,j] = ginv[i,j].subs( g0'+str(p)+'_'+str(2*k+1)+' == f_g0'+str(p)+'_'+str(2*k+1)+')')
    
if show==True:
    display('g = ',g)
    display(r'g^{-1} = ',ginv)
                


# Una vez disponemos tanto de la métrica como de su inversa expresados únicamente en términos de las funciones $g_{\mu\nu}^{(k)}$, podemos calcular la conexión afín mediante su relación con la métrica:
# 
# $$\Gamma^{\mu}_{\nu\lambda} = \frac{1}{2} g^{\mu\rho}\left\{ \frac{\partial g_{\nu\rho}}{\partial x^{\lambda}} + \frac{\partial g_{\lambda\rho}}{\partial x^{\nu}} - \frac{\partial g_{\nu\lambda}}{\partial x^{\rho}}  \right\}$$
# 
# Dado que los elementos de la conexión afín los necesitamos para computar la aceleración de una partícula
# 
# $$\frac{d^2 x^i}{dt^2} = -\Gamma^i_{00} - 2\Gamma^i_{0j} \frac{dx^j}{dt} - \Gamma^i_{jk} \frac{dx^j}{dt}\frac{dx^k}{dt} + \left[ \Gamma^0_{00} + 2\Gamma^0_{0j}\frac{dx^j}{dt} + \Gamma^0_{jk} \frac{dx^j}{dt}\frac{dx^k}{dt} \right]\frac{dx^i}{dt}$$
# 
# truncamos cada término al orden de potencias de $v$ necesario para que la aceleración quede a orden $n$.

# In[10]:


# Definimos la conexión afín

Gamma0 = Matrix(4,aux)
Gamma1 = Matrix(4,aux)
Gamma2 = Matrix(4,aux)
Gamma3 = Matrix(4,aux)

for mu in range(4):
    for nu in range(4):
        for lbda in range(4):
            for rho in range(4):
                exec('Gamma'+str(mu)+'[nu,lbda] += 1/2 * ginv[mu,rho] * ( diff_v(g[rho,nu], lbda) + diff_v(g[rho,lbda], nu) - diff_v(g[nu,lbda], rho ) ) ')                
            
Gamma0 = Gamma0(aux=0)
Gamma1 = Gamma1(aux=0)
Gamma2 = Gamma2(aux=0)
Gamma3 = Gamma3(aux=0)

# Truncamos cada término en el necesario para obtener orden n en v

#gamma000 a orden n-1
Gamma0[0,0] = truncate(Gamma0[0,0], v, n-1)

#gammai00 a orden n
for i in range(1,4):
    exec('Gamma'+str(i)+'[0,0] = truncate(Gamma'+str(i)+'[0,0], v, n)')

#gamma00j a orden n-2
for j in range(1,4):
    exec('Gamma0[0,j] = truncate(Gamma0[0,j], v, n-2)')
    exec('Gamma0[j,0] = truncate(Gamma0[j,0], v, n-2)')
    
#gamma0jk a orden n-3 
for j in range(1,4):
    for k in range(1,4):
        exec('Gamma0['+str(j)+','+str(k)+'] = truncate(Gamma0['+str(j)+','+str(k)+'], v, n-3)')
        
#gammai0j a orden n-1
for i in range(1,4):
    for j in range(1,4):
        exec('Gamma'+str(i)+'[0,j] = truncate(Gamma'+str(i)+'[0,j], v, n-1)')
        exec('Gamma'+str(i)+'[j,0] = truncate(Gamma'+str(i)+'[j,0], v, n-1)')
        
#gammaijk a orden n-2
for i in range(1,4):
    for j in range(1,4):
        for k in range(1,4):
            exec('Gamma'+str(i)+'[j,k] = truncate(Gamma'+str(i)+'[j,k], v, n-2)')

            
if(show==True):
    for i in range(4):
        for j in range(4):
            for k in range(j,4):
                exec('display("Gamma("+str(i)+","+str(j)+","+str(k)+") = ")')
                exec('display(Gamma'+str(i)+'['+str(j)+','+str(k)+'])')


# Una vez tenemos la métrica y la conexión afín, podemos calcular el tensor de Ricci
# 
# $$R_{\mu\kappa}=\frac{\partial \Gamma^\lambda_{\mu\lambda}}{\partial x^\kappa}-\frac{\partial\Gamma^\lambda_{\mu\kappa}}{\partial x^\lambda}+\Gamma^\eta_{\mu\lambda} \Gamma^\lambda_{\kappa\eta} - \Gamma^\eta_{\mu\kappa} \Gamma^\lambda_{\eta\lambda}$$
# 
# Dado que, para los coeficientes de la conexión afín de los que disponemos (truncados a los órdenes de magnitud indicados anteriormente), solo podemos calcular completamente los términos del tensor de Ricci
# 
# $$R_{00} = R_{00}^{(2)} + R_{00}^{(4)} + ... + R_{00}^{(n)}$$
# 
# $$R_{ij} = R_{ij}^{(2)} + R_{ij}^{(4)} + ... + R_{ij}^{(n-2)}$$
# 
# $$R_{0i} = R_{0i}^{(3)} + R_{0i}^{(5)} + ... + R_{0i}^{(n-1)}$$
# 
# truncamos cada término al orden correspondiente.

# In[11]:


Ricci = Matrix(4, aux)

for mu in range(4):
    for kappa in range(mu,4):
        for lbda in range(4):
            exec('Ricci[mu,kappa] += diff_v( Gamma'+str(lbda)+'[mu,lbda], kappa ) - diff_v( Gamma'+str(lbda)+'[mu,kappa], lbda )')
            for eta in range(4):
                exec('Ricci[mu,kappa] += Gamma'+str(eta)+'[mu,lbda]*Gamma'+str(lbda)+'[kappa,eta] - Gamma'+str(eta)+'[mu,kappa]*Gamma'+str(lbda)+'[eta,lbda]')
                
Ricci = Ricci(aux=0)


# R00 a orden n
Ricci[0,0] = truncate(Ricci[0,0], v, n)
        
# Rij a orden n-2
for i in range(1,4):
    for j in range(i,4):
        polinomio=Ricci[i,j]
        Ricci[i,j] = Ricci[j,i] = truncate(Ricci[i,j], v, n-2)
        
# R0i a orden n-1
for i in range(1,4):
    Ricci[0,i] = Ricci[i,0] = truncate(Ricci[0,i], v, n-1)

if(show==True):
    for i in range(4):
        for j in range(i,4):
            display('Ricci('+str(i)+','+str(j)+')')
            display(Ricci[i,j].full_simplify())
        


# Para obtener las ecuaciones del movimiento calculamos las aceleraciones $d^2x^i/dt^2$ a través de la conexión afín como
# 
# $$\frac{d^2 x^i}{dt^2} = -\Gamma^i_{00} - 2\Gamma^i_{0j} \frac{dx^j}{dt} - \Gamma^i_{jk} \frac{dx^j}{dt}\frac{dx^k}{dt} + \left[ \Gamma^0_{00} + 2\Gamma^0_{0j}\frac{dx^j}{dt} + \Gamma^0_{jk} \frac{dx^j}{dt}\frac{dx^k}{dt} \right]\frac{dx^i}{dt}$$
# 

# In[12]:


def runge_kutta(h, vec, eq1, eq2, eq3):
    vec_new = [0,0,0,0,0,0,0]
    
    KL0 = [[0,0] for i in range(3)]  # K^i_j = KLj[i][0] ; L^i_j = KLj[i][1]
    KL1 = [[0,0] for i in range(3)]
    KL2 = [[0,0] for i in range(3)]
    KL3 = [[0,0] for i in range(3)]
    
    for i in range(1,4):
        exec('KL0[i-1][0] = eq'+str(i)+'(x1=vec[1], x2=vec[2], x3=vec[3], v1=vec[4], v2=vec[5], v3=vec[6]).n()')
        KL0[i-1][1] = vec[i+3]
    
    for i in range(1,4):
        exec('KL1[i-1][0] = eq'+str(i)+'(x1=vec[1]+0.5*KL0[0][1]*h, x2=vec[2]+0.5*KL0[1][1]*h, x3=vec[3]+0.5*KL0[2][1]*h, v1=vec[4]+0.5*KL0[0][0]*h, v2=vec[5]+0.5*KL0[1][0]*h, v3=vec[6]+0.5*KL0[2][0]*h).n()')
        KL1[i-1][1] = vec[i+3]+0.5*KL0[i-1][0]*h
    
    for i in range(1,4):
        exec('KL2[i-1][0] = eq'+str(i)+'(x1=vec[1]+0.5*KL1[0][1]*h, x2=vec[2]+0.5*KL1[1][1]*h, x3=vec[3]+0.5*KL1[2][1]*h, v1=vec[4]+0.5*KL1[0][0]*h, v2=vec[5]+0.5*KL1[1][0]*h, v3=vec[6]+0.5*KL1[2][0]*h).n()')
        KL1[i-1][1] = vec[i+3]+0.5*KL1[i-1][0]*h

    for i in range(1,4):
        exec('KL3[i-1][0] = eq'+str(i)+'(x1=vec[1]+KL2[0][1]*h, x2=vec[2]+KL2[1][1]*h, x3=vec[3]+KL2[2][1]*h, v1=vec[4]+KL2[0][0]*h, v2=vec[5]+KL2[1][0]*h, v3=vec[6]+KL2[2][0]*h).n()')
        KL1[i-1][1] = vec[i+3]+KL2[i-1][0]*h
        
    vec_new[0] = vec[0] + h
    vec_new[1] = vec[1] + 1/6.0 * (KL0[0][1] + 2*KL1[0][1] + 2*KL2[0][1] + KL3[0][1]) * h
    vec_new[2] = vec[2] + 1/6.0 * (KL0[1][1] + 2*KL1[1][1] + 2*KL2[1][1] + KL3[1][1]) * h
    vec_new[3] = vec[3] + 1/6.0 * (KL0[2][1] + 2*KL1[2][1] + 2*KL2[2][1] + KL3[2][1]) * h
    vec_new[4] = vec[4] + 1/6.0 * (KL0[0][0] + 2*KL1[0][0] + 2*KL2[0][0] + KL3[0][0]) * h
    vec_new[5] = vec[5] + 1/6.0 * (KL0[1][0] + 2*KL1[1][0] + 2*KL2[1][0] + KL3[1][0]) * h
    vec_new[6] = vec[6] + 1/6.0 * (KL0[2][0] + 2*KL1[2][0] + 2*KL2[2][0] + KL3[2][0]) * h
    
    return vec_new


# In[13]:


#Ecuaciones del movimiento

v1 = var('v1', latex_name='\\frac{dr}{dt}')
v2 = var('v2', latex_name='\\frac{d\\theta}{dt}')
v3 = var('v3', latex_name='\\frac{d\\phi}{dt}')

eq_a1 = eq_a2 = eq_a3 = 0
eq_v1 = v1
eq_v2 = v2
eq_v3 = v3

for i in range(1,4):
    exec('eq_a'+str(i)+' += -Gamma'+str(i)+'[0,0] + Gamma0[0,0]*v'+str(i))
    for j in range(1,4):
        exec('eq_a'+str(i)+' += -2*Gamma'+str(i)+'[0,j]*v'+str(j)+' + 2*Gamma0[0,j]*v'+str(j)+'*v'+str(i))
        for k in range(1,4):
            exec('eq_a'+str(i)+' += -Gamma'+str(i)+'[j,k]*v'+str(j)+'*v'+str(k)+' + Gamma0[j,k]*v'+str(j)+'*v'+str(k)+'*v'+str(i))
    exec('eq_a'+str(i)+' = eq_a'+str(i)+'(v=1)')

#eq_a1 += 2*G*M*v1**2/c**2/x1**2

display('a1 = ',eq_a1)
display('a2 = ',eq_a2)
display('a3 = ',eq_a3)


# In[14]:


#Constantes
const_G=1     #6.67*10**(-11)
const_M=1     #2*10**(30)
const_c=1     #3*10**8

for i in range(1,4):
    exec('eq_a'+str(i)+' = eq_a'+str(i)+'(G=const_G, M=const_M, c=const_c)')

#Condiciones iniciales [t,x1,x2,x3,v1,v2,v3], tiempo de simulación y paso temporal
ics = [0, 100, acos(-1)/2+0.2, 0, 0, 0.005, 0.000000005]
tfin = 1000
deltat = 1

sol = [[0,0,0,0,0,0,0] for i in range(int(tfin/deltat))]
sol[0] = ics

for i in range(1,len(sol)):
    if i%100==0:
        print(i)
    sol[i] = runge_kutta(deltat, sol[i-1], eq_a1, eq_a2, eq_a3)


# In[15]:


# XYZ

S = [[i,j,k] for a,i,j,k,b,c,d in sol]
S_cart = [[0,0,0] for i in range(len(S))]

for i in range(len(S)):
    S_cart[i][0] = S[i][0]*sin(S[i][1])*cos(S[i][2])
    S_cart[i][1] = S[i][0]*sin(S[i][1])*sin(S[i][2])
    S_cart[i][2] = S[i][0]*cos(S[i][1])

T = point3d(S_cart)

S = Sphere(2*const_G*const_M/const_c**2, color='orange')

(S+T).show()

