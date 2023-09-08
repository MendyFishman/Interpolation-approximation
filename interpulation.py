import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import Symbol
from sympy import simplify, expand
import sympy as sp
from sympy import sin, cos, pi
from fractions import Fraction


#Mendy Fishman 
#Lev Velitski 


def compute_points(a,b,n):
    points = np.zeros(n+1)
    for k in range(0,n+1):
        points[k] = a + ((b-a)/n)*k
    return points

def chebichev_points(a,b,n):
    #compute points in [-1,1]:
    zeta_array = [0]*n
    for k in range(1,n+1):
        expression_in_zeta_cos = ((2*k-1)/(2*n))*pi
        zeta_array[k-1] = cos(expression_in_zeta_cos)
        
    #transform from [-1,1] -> [a,b] by z = m1*x+n1:
    #z[-1] = -m1 + n1 = a
    #z[1] = m1 + n1 = b
    #So, we get:  n1 = (a+b)/2 ,  m1 = (b-a)/2 
    m1 = (b-a)/2 
    n1 = (a+b)/2
    for k in range(0,n):
        zeta_array[k] = m1*zeta_array[k] + n1
    for k in range(0,n):
        zeta_array[k] = float(zeta_array[k])
    return zeta_array
    
    
def compute_elementary_polynomial(points):
    n = len(points)
    x = Symbol('x')
    l = [0]*n
    t = 0
    for k in range(0,n):
        l_k = 1
        for j in range(0,n):
            if k!=j:
                t = (x-points[j])/(points[k]-points[j])
                l_k = t*l_k
        l[k] = expand(l_k)       
    return l
    
def compute_approximation(l,f, points):
        x = Symbol('x')
        n = len(points)
        L_n = 0
        for i in range(0,n):
            L_n += f.subs(x,points[i])*l[i]
        return L_n
    
def plot_functions(f,Ln,a,b,points,L_n_cheb,Chebichev_points):
    x = Symbol('x')
    n = len(points)
    n_20 = 20*n
    m = np.linspace(a, b, num = n_20)
    f_array = np.zeros(n_20)
    Ln_array = np.zeros(n_20)
    Ln_Cheb_array = np.zeros(n_20)
    for i in range(0,n_20):
        f_array[i] = f.subs(x,m[i])
        Ln_array[i] = Ln.subs(x,m[i])
        Ln_Cheb_array[i] = L_n_cheb.subs(x,m[i]) 
    # Plotting both the curves simultaneously
    plt.plot(m, f_array, color='r', label='f = ' + str(f))
    plt.plot(m, Ln_array, color='g', label='L'+str(n-1))
    plt.plot(m, Ln_Cheb_array, color='b', label='L_Cheb'+str(n-1))
    
    # Naming the x-axis, y-axis and the whole graph
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Interpulation Approximation (by Lagrange)" + " n = " + str(n-1))
    plt.legend()
    plt.show()

def Get_Approximation(f,a,b,n):
    x = Symbol('x')
    
    #compute the division points:
    points_equal_distance = compute_points(a,b,n)

    
    #compute Chebichev points:
    Chebichev_points = chebichev_points(a,b,n+1)
    Chebichev_points.sort()
    
    #print the points arrays:
    #equal distance points:
    print("a = " + str(a) + ", b = " + str(b) + ", n = " + str(n))
    print("") 
    print("The points of equal distance division: ",points_equal_distance)
    print("") 
    #Chebichev points:
    print("The points by Chebichev: ",Chebichev_points)
    print("")    
    
    #compute the lagrange elementary polynomial:
    l = compute_elementary_polynomial(points_equal_distance)
    print("The elementary polynomial: ")
    print(l)
    print("")
    
    #compute the lagrange elementary polynomial by Chebichev points:
    l_cheb = compute_elementary_polynomial(Chebichev_points)
    print("The elementary polynomial (by Chebicev points): ")
    print(l_cheb)
    print("")
    
    #compute the approximation by equal distance points:
    L_n = compute_approximation(l,f, points_equal_distance)
    print("The interpulation approximation by equal diance points is: ")
    print("L_" + str(n) + " = ", L_n)
    
    #compute the approximation by chebichev points:
    L_n_cheb = compute_approximation(l_cheb,f, Chebichev_points)
    print("The interpulation approximation by Chebichev points is: ")
    print("L_" + str(n) + " = ", L_n_cheb)
        
    #plot the function and the approximation:
    plot_functions(f,L_n,a,b,points_equal_distance,L_n_cheb,Chebichev_points) 




def main():
    x = Symbol('x')
    
    #The inputs:
    f =  x**0.5/(x+1)         
    a = 0
    b = 1
    n = 2
    
    #activate the function
    Get_Approximation(f,a,b,n)

main()


