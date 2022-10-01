from sympy import *
import copy

#Definiciones
M,t,r,theta,phi ,e = symbols('M t r theta phi e')
x,y,z = symbols('x y z')
v=Function('v')(t)
#define the coordinates
coord = [t,x,y,z ]
coord_nombres = ["t","r","\theta","\phi"]
I2=eye(4)
I3=eye(4)

I3[0,0]=-1

#Orden hasta la que realizar la expansión en serie en términos de velocidad
orden=3

orden+=1

if orden %2 ==0:
    orden_par=orden
else:
    orden_par=orden+1
    
G=zeros(4)
def delta_dirac(i,j):
    if i==j:
        return 1
    else:
        return 0
    
m2=open("metrica2.txt","r")
M2=zeros(4)
M2[0,0]=t**2
M2[1,1]=2*x
M2[2,2]=3*y
M2[3,3]=4*z

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r2=M2[{i},{j}]*v**2")
        
M3=zeros(4)
M3[0,1]=M3[1,0]=2*x
M3[0,2]=M3[2,0]=3*y
M3[0,3]=M3[3,0]=4*z

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r3=M3[{i},{j}]*v**3")
        
M4=zeros(4)
M4[0,0]=t**5
M4[1,1]=2*x**2
M4[2,2]=3*y**2
M4[3,3]=4*z**2

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r4=M4[{i},{j}]*v**4")
        
M5=zeros(4)
M5[0,1]=M5[1,0]=2*x**3
M5[0,2]=M5[2,0]=3*y**3
M5[0,3]=M5[3,0]=4*z**3

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r5=M5[{i},{j}]*v**5")
        

#Definición de los coeficientes de la métrica como funciones
m2=open("metrica2.txt","r")
M2=zeros(4)
M2[0,0]=g002f=Function('g_00^2')(x,y,z,t)
M2[1,1]=Function('g_11^2')(x,y,z,t)
M2[2,2]=Function('g_22^2')(x,y,z,t)
M2[3,3]=Function('g_33^2')(x,y,z,t)

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r2=M2[{i},{j}]*v**2")
        exec(f"g{i}{j}r0=I2[{i},{j}]")
        exec(f"g{i}{j}r1=0")
        
        
M3=zeros(4)
M3[0,1]=M3[1,0]=Function('g_01^3')(x,y,z,t)
M3[0,2]=M3[2,0]=Function('g_02^3')(x,y,z,t)
M3[0,3]=M3[3,0]=Function('g_03^3')(x,y,z,t)

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r3=M3[{i},{j}]*v**3")
        
M4=zeros(4)
M4[0,0]=Function('g_00^4')(x,y,z,t)
M4[1,1]=Function('g_11^4')(x,y,z,t)
M4[2,2]=Function('g_22^4')(x,y,z,t)
M4[3,3]=Function('g_33^4')(x,y,z,t)

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r4=M4[{i},{j}]*v**4")
        
M5=zeros(4)
M5[0,1]=M5[1,0]=Function('g_00^5')(x,y,z,t)
M5[0,2]=M5[2,0]=Function('g_00^5')(x,y,z,t)
M5[0,3]=M5[3,0]=Function('g_00^5')(x,y,z,t)

for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}r5=M5[{i},{j}]*v**5")
            
        
#Definimos como símbolos los elementos de la métrica y sus coeficientes de la expansión en serie
for i in range(4):
    for j in range(4):
        exec(f"g{i}{j}=symbols('g_{i}{j}')")    
        #exec(f"g{i}{j}r=symbols('g_{i}{j}')")    

        #exec(f"g{i}{j}r=Function('g_{i}{j}r')(t,x,y,z,v)")    
        exec(f"g{i}{j}={I3[i,j]}")        
        for k in range(orden):
            exec(f"g{i}{j}{k} = symbols('g^{k}_{i}{j}')")

            #exec(f"g{i}{j}{k}r = Function('g^{k}_{i}{j}r')(t,x,y,z,v)")

#Escribimos explícitamente la expansión en serie
for i in range(4):
    for j in range(i,4):       
        for k in range(1,orden):
            exec(f"g{i}{j}=g{i}{j}+e**{k}*g{i}{j}{k}")
            #exec(f"g{i}{j}r=g{i}{j}r+e**{k}*g{i}{j}r{k}")
                       
            
            
            
#Órdenes que se anulan por la naturaleza del problema
for i in range(1,4):
    for k in range(2,orden_par-1,2):
        exec(f"g0{i}=g0{i}.subs(g0{i}{k},0)")           
    exec(f"g0{i}=g0{i}.subs(g0{i}1,0)")
for i in range(1,4):
    for j in range(i,4):
        for k in range(1,orden,2):
            exec(f"g{i}{j}=g{i}{j}.subs(g{i}{j}{k},0)")  
            
            exec(f"g00=g00.subs(g00{k},0)")




for i in range(4):
    for j in range(0,i):
        exec(f"g{i}{j}=g{j}{i}")

        
#Generamos la matriz de la métrica        
G=zeros(4)
for i in range(4):
    for j in range(i,4):
        exec(f"G[{i},{j}]=g{i}{j}")
for i in range(4):
    for j in range(0,i):
        exec(f"G[{i},{j}]=g{j}{i}")
        
        
#Definimos como símbolos los elementos de la inversa de la métrica y sus coeficientes de la expansión en serie
for i in range(4):
    for j in range(4):
        exec(f"gg{i}{j}=symbols('g^{i}{j}')")    
        exec(f"gg{i}{j}r=Function('g^{i}{j}r')(t,x,y,z,v)")    
        exec(f"gg{i}{j}={I3[i,j]}")        
        for k in range(orden):
            exec(f"gg{i}{j}{k} = symbols('g^{i}{j}{k}')")
            exec(f"gg{i}{j}{k}r = Function('g^{i}{j}{k}r')(t,x,y,z,v)")

#Escribimos explícitamente la expansión en serie
for i in range(4):
    for j in range(i,4):       
        for k in range(1,orden):
            exec(f"gg{i}{j}=gg{i}{j}+e**{k}*gg{i}{j}{k}")
            exec(f"gg{i}{j}r=gg{i}{j}r+e**{k}*gg{i}{j}{k}r")
            
            
#Órdenes que se anulan por su naturaleza
for i in range(1,4):
    for k in range(2,orden-1,2):
        exec(f"gg0{i}=gg0{i}.subs(gg0{i}{k},0)")           
    exec(f"gg0{i}=gg0{i}.subs(gg0{i}1,0)")
for i in range(1,4):
    for j in range(i,4):
        for k in range(1,orden,2):
            exec(f"gg{i}{j}=gg{i}{j}.subs(gg{i}{j}{k},0)")           
for k in range(1,orden_par-1,2):
    exec(f"gg00=gg00.subs(gg00{k},0)")    
    




#Definición de ciertas funciones útiles para usos posteriores
def cortar_orden(expr,orden=4):
    coef=zeros(1,orden+1)
    expresion=0
    for i in range(orden+1):
        coef[i]=expr.expand().coeff(e,i)
        expresion+=coef[i]*e**i
    return expresion

def dejar_orden(expr,orden=4):
    return expr.expand().coeff(e,orden)

def cortar_orden_matriz(matr,orden=4,dim=4):
    m=zeros(dim)
    for i in range(dim):
        for j in range(dim):
            m[i,j]=cortar_orden(matr[i,j],orden)
    return m

def dejar_orden_matriz(matr,orden=4,dim=4):
    return cortar_orden_matriz(matr,orden+1,dim)-cortar_orden_matriz(matr,orden-1,dim)

def sustituir_matriz(matr,simbolo,valor,dim=4):
    m=copy(matr)
    for i in range(dim):
        for j in range(dim):
            m[i,j]=m[i,j].subs(simbolo,valor)
            
    return m 

#Definimos la matriz inversa
Ginv=zeros(4)
for ii in range(4):
    for jj in range(ii,4):
        exec(f"Ginv[{ii},{jj}]=gg{ii}{jj}")
for ii in range(4):
    for jj in range(0,ii):
        exec(f"Ginv[{ii},{jj}]=gg{jj}{ii}")
        
##Con las ecuaciones obtenidas en I1, resolvemos el sistema y sustituimos los terminos de Ginv en función de los de G, obteniendo Ginv2
Ginv2=zeros(4)
I1=Ginv*G
I1=cortar_orden_matriz(I1.expand(),orden-1)-I2
espacio=" /  "
for i in range(4):
    for j in range(i,4):
        for k in range(1,orden):
            for k2 in range(1,orden):
               
                exec(f"cond=(len(solve((I1[{i},{j}]).coeff(e,{k2}),gg{i}{j}{k}))>0)")
               
                if cond:
                    exec(f"h{i}{j}{k}=solve((I1[{i},{j}]).coeff(e,{k2}),gg{i}{j}{k})[0]")  ##Resolvemos las ecuaciones                
                    for iii in range(4):
                        for jjj in range(iii,4):
                            exec(f"gg{iii}{jjj}=gg{iii}{jjj}.subs(gg{i}{j}{k},h{i}{j}{k})") ##Sustituimos la solución 
                                                                                            ##en todos los términos de gg
                                
                            exec(f"I1[{iii},{jjj}]=I1[{iii},{jjj}].subs(gg{i}{j}{k},h{i}{j}{k})") ##Sustituimos la solución en 
                                                                                                  ##todos los términos de I1
             

            
##Una vez resuelto el sistema, solo queda igualar término a término la matriz resultante
for ii in range(4):
    for jj in range(ii,4):
        exec(f"Ginv2[{ii},{jj}]=gg{ii}{jj}")
for ii in range(4):
    for jj in range(0,ii):
        exec(f"Ginv2[{ii},{jj}]=gg{jj}{ii}")
       
#Sustitución de los símbolos por funciones

Ginv3=copy.deepcopy(Ginv2)               
for i in range(4):
    for j in range(4):
        for k in range(2,orden):
            exec(f"Ginv3=Ginv3.subs(g{i}{j}{k},g{i}{j}r{k})")
for ii in range(4):
    for jj in range(ii,4):
        exec(f"Ginv[{ii},{jj}]=gg{ii}{jj}")
for ii in range(4):
    for jj in range(0,ii):
        exec(f"Ginv[{ii},{jj}]=gg{jj}{ii}")
        
#Definición de los símbolos de las derivadas parciales de las componentes de la métrica
for i in range(4):
    for j in range(4):
        exec(f"d{k}g{i}{j}=0")
        for k in range(4):
            exec(f"d{k}g{i}{j}=symbols('dg_{coord[k]}{i}{j}')")
            exec(f"d{k}g{i}{j}=0")
            for h in range(orden):
                exec(f"d{k}g{i}{j}{h}=symbols('dg_{coord[k]}{i}{j}^{h}')")
                if k==0:
                    exec(f"d{k}g{i}{j}=d{k}g{i}{j}+e**{h+1}*d{k}g{i}{j}{h}")
                else:
                    exec(f"d{k}g{i}{j}=d{k}g{i}{j}+e**{h}*d{k}g{i}{j}{h}")
         
        
for h in range(4):
    exec(f"d{h}G=zeros(4)") 
    for i in range(4):
        for j in range(i,4):
            exec(f"d{h}G[{i},{j}]=d{h}g{i}{j}")
    for i in range(4):
        for j in range(0,i):
            exec(f"d{h}G[{i},{j}]=d{h}g{j}{i}")
     
#Órdenes que se anulan por la naturaleza del problema    
for h in range(4):         
    for i in range(1,4):
        for k in range(0,orden,2):
            exec(f"d{h}g0{i}=d{h}g0{i}.subs( d{h}g0{i}{k},0)")           
        exec(f"d{h}g0{i}=d{h}g0{i}.subs(d{h}g0{i}1,0)")
        exec(f"d{h}G[0,{i}]=d{h}G[0,{i}].subs(d{h}g0{i}1,0)")
        exec(f"d{h}G[0,{i}]=d{h}G[0,{i}].subs(d{h}g0{i}0,0)")
    for i in range(1,4):
        for j in range(4):
            for k in range(1,orden,2):
                exec(f"d{h}g{i}{j}=d{h}g{i}{j}.subs(d{h}g{i}{j}{k},0)")   
            exec(f"d{h}g{i}{j}=d{h}g{i}{j}.subs(d{h}g{i}{j}0,0)")     
            exec(f"d{h}g00=d{h}g00.subs(d{h}g00{k},0)")
    for i in range(4):
        for j in range(i,4):
            exec(f"d{h}G[{i},{j}]=d{h}g{i}{j}")
    for i in range(4):
        for j in range(0,i):
            exec(f"d{h}G[{i},{j}]=d{h}g{j}{i}")
        
#Definición de las funciones de las derivadas parciales
for i in range(4):
    for j in range(4):
        for k in range(4):
            for h in range(0,orden):
                exec(f"d{k}g{i}{j}r{h}=diff(g{i}{j}r{h},coord[{k}])")



#Definición y obtención de los coeficientes de la conexión afín
for i in range(4):
    for j in range(4):
        for k in range(j,4): 
            exec(f"T{i}{j}{k}=symbols('\Gamma^{i}_{j}{k}')")
            exec(f"T{i}{j}{k}=0")
            for h in range(4):
                exec(f"T{i}{j}{k}+=0.5*Ginv2[i,h]*(d{k}G[{h},{j}]+d{j}G[{h},{k}]-d{h}g{j}{k})")
            exec(f"T{i}{j}{k}=cortar_orden(T{i}{j}{k},orden-2)")

            exec(f"T{i}{k}{j}=T{i}{j}{k}")
            
for i in range(4):
    for j in range(4):
        for k in range(4):
            for h in range(orden):
                exec(f"T{i}{j}{k}{h}=symbols('\Gamma^{i}_{j}{k}^{h}')")
                exec(f"T{i}{j}{k}{h}=T{i}{j}{k}.expand().coeff(e,{h})")
                
#Sustitución de los símbolos por funciones en las componentes de la conexión afín
for ii in range(4):
    for jj in range(4):
        for kk in range(4):
            exec(f"T{ii}{jj}{kk}r=Function('T{ii}{jj}{kk}r')")
            exec(f"T{ii}{jj}{kk}r=T{ii}{jj}{kk}")

                        
for ii in range(4):
    for jj in range(4): 
        for kk in range(4):
            for i in range(4):
                for j in range(4):
                    for k2 in range(orden):
                        exec(f"T{ii}{jj}{kk}r=T{ii}{jj}{kk}r.subs(g{i}{j}{k2},g{i}{j}r{k2})")
                        for k in range(4):
                            exec(f"T{ii}{jj}{kk}r=T{ii}{jj}{kk}r.subs(d{k}g{i}{j}{k2},d{k}g{i}{j}r{k2})") 

                        
#Sustitución de los símbolos por funciones en las componentes de la conexión afín
for ii in range(4):
    for jj in range(4):
        for kk in range(4):
            for h in range(orden):
                exec(f"T{ii}{jj}{kk}r{h}=Function('T{ii}{jj}{kk}r{h}')")
                exec(f"T{ii}{jj}{kk}r{h}=T{ii}{jj}{kk}r.coeff(e,{h})")
#Definición de las funciones "derivadas parciales de las componentes de la conexión afín"
for i in range(4):
    for j in range(4):
        for k in range(4):
            for h in range(4):
                for m in range(orden):
                    exec(f"d{k}T{i}{j}{h}r{m}=diff(T{i}{j}{h}r.coeff(e,{m}),coord[{k}])")                 
for f in range(4):
    for i in range(4):
        for j in range(4):
            for k in range(4):
                exec(f"d{k}T{f}{i}{j}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}')")
                exec(f"d{k}T{f}{i}{j}=0")
                for h in range(orden):
                    exec(f"d{k}T{f}{i}{j}{h}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}^{h}')")
                    exec(f"cond=(T{f}{i}{j}.coeff(e,{h})==0)")
                    if cond:
                         exec(f"d{k}T{f}{i}{j}{h}=0")
                    if k==0:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h+1}*d{k}T{f}{i}{j}{h}")
                    else:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h}*d{k}T{f}{i}{j}{h}")
                        
#Definición de las funciones "derivadas parciales de las componentes de la conexión afín"
for i in range(4):
    for j in range(4):
        for k in range(4):
            for h in range(4):
                for m in range(orden):
                    exec(f"d{k}T{i}{j}{h}r{m}=diff(T{i}{j}{h}r.coeff(e,{m}),coord[{k}])")                 
for f in range(4):
    for i in range(4):
        for j in range(4):
            for k in range(4):
                exec(f"d{k}T{f}{i}{j}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}')")
                exec(f"d{k}T{f}{i}{j}=0")
                for h in range(orden):
                    exec(f"d{k}T{f}{i}{j}{h}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}^{h}')")
                    exec(f"cond=(T{f}{i}{j}.coeff(e,{h})==0)")
                    if cond:
                         exec(f"d{k}T{f}{i}{j}{h}=0")
                    if k==0:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h+1}*d{k}T{f}{i}{j}{h}")
                    else:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h}*d{k}T{f}{i}{j}{h}")
                        
#Definición de las funciones "derivadas parciales de las componentes de la conexión afín"
for i in range(4):
    for j in range(4):
        for k in range(4):
            for h in range(4):
                for m in range(orden):
                    exec(f"d{k}T{i}{j}{h}r{m}=diff(T{i}{j}{h}r.coeff(e,{m}),coord[{k}])")                 
for f in range(4):
    for i in range(4):
        for j in range(4):
            for k in range(4):
                exec(f"d{k}T{f}{i}{j}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}')")
                exec(f"d{k}T{f}{i}{j}=0")
                for h in range(orden):
                    exec(f"d{k}T{f}{i}{j}{h}=symbols('d\Gamma_{coord[k]}^{f}_{i}{j}^{h}')")
                    exec(f"cond=(T{f}{i}{j}.coeff(e,{h})==0)")
                    if cond:
                         exec(f"d{k}T{f}{i}{j}{h}=0")
                    if k==0:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h+1}*d{k}T{f}{i}{j}{h}")
                    else:
                        exec(f"d{k}T{f}{i}{j}=d{k}T{f}{i}{j}+e**{h}*d{k}T{f}{i}{j}{h}")
                        
                        
#Sustitución de los símbolos por funciones del tensor de Ricci
for ii in range(4):
    for jj in range(4):
        exec(f"R{ii}{jj}r=Function('R{ii}{jj}r')")
        exec(f"R{ii}{jj}r=R{ii}{jj}")
for ii in range(4):
    for jj in range(4):            
        for i in range(4):
            for j in range(4):
                for h in range(orden):
                    exec(f"R{ii}{jj}r=R{ii}{jj}r.subs(g{i}{j}{h},g{i}{j}r{h})")

                        
for ii in range(4):
    for jj in range(4):                      
        for i in range(4):
            for j in range(4):
                for h in range(4):        
                    for k in range(4):
                        for kk in range(orden):
                            exec(f"R{ii}{jj}r=R{ii}{jj}r.subs(d{k}g{i}{j}{kk},d{k}g{i}{j}r{kk})") 
                            exec(f"R{ii}{jj}r=R{ii}{jj}r.subs(d{k}T{i}{j}{h}{kk},d{k}T{i}{j}{h}r{kk})") 
