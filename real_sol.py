import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from matplotlib import pyplot as plt
def rwf(a,b,i):
    while a<b:
        yield a
        a+=i
max_rect=0
e=2.718
gamma=[]
re=[]
t_f=[]
t_r=[]

for a in range(1,100):

    def k1(x):
        return 2*(e**(0.0001*a*x))+a
    def k2(x):
        return (2*(e**(-0.0001*a*(x-800))))
    max_rect=0
    max_gamma=0
    Tf_max=0
    Tr_max=0

    for g in rwf(0,2,0.0001):
        
        def func1(tf1):
            res1,err1= quad(k1,tf1,500)
            res2,err2= quad(k2,273,tf1)
            res=(g*res1)-res2
            return res

        # Use 1.0 as the initial guess.  Note that a bad initial guess
        # might generate a warning and return the degenerate solution at x=0.
        tf = fsolve(func1, 273)

        def func2(tr1):
            res1,err1=quad(k2,tr1,500)
            res2,err2=quad(k1,273,tr1)
            res=res1-(g*res2)
            return res
        

        tr=fsolve(func2,273)

        r1,e1=quad(k1,tf,500)
        r2,e2=quad(k1,273,tr)
        rect=r1/r2

        if rect>max_rect:
            max_rect =rect
            max_gamma=g
            Tf_max=tf
            Tr_max=tr

        #print (g,tf,tr,rect)
    gamma.append(max_gamma)
    re.append(max_rect)
    t_f.append(Tf_max)
    t_r.append(Tr_max)
        
    #print("\n\nMaximum rectification ratio=",max_rect,"\nMaximum gamma=",max_gamma,"\nTf_max=",Tf_max,"\nTr_max=",Tr_max)



#plt.plot(gamma,re)
plt.plot(t_f,t_r)

plt.xlabel("T_f")
plt.ylabel("T_r")

plt.show()
