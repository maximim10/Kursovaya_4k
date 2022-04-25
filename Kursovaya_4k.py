import math
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool

L=1; n=20; h=L/n; b=h/1000; gamma=1/2; E1=2; E2=E1/10; points=100000; scale=n/3
h0 = 1/(gamma/(E1)+(1-gamma)/(E2)); c2 = 1/2-h0*(gamma**2/E1/2+(1-gamma)*gamma/E1+(1-gamma)**2/E2/2)
k1=0
def Find_force(t):
    return 2*b**2/((t-L/2)**2+b**2)
def Find_G_part(t):
    G=b**2*math.log(4*b**2+(L-2*t)**2)-b*(L-2*t)*math.atan((L-2*t)/(2*b))
    return G
def Find_G(t, k1):
    G=Find_G_part(t)+k1*t
    return G
def Find_k1():
    global k1
    summ=0
    for k in range (0,n):
        summ+=(Find_G_part(k*h+gamma*h)-Find_G_part(k*h))/E1
        summ+=(Find_G_part(k*h+h)-Find_G_part(k*h+gamma*h))/E2
    k1=-summ/Sum_near_k1()
    return k1
def Sum_near_k1():
    summ=0
    for k in range (0,n):
        summ+=gamma*h/E1
        summ+=h*(1-gamma)/E2
    return summ
def Get_k1():
    return k1

def displacement_exact(t):
    d=0
    for k in range (0,n):
        if (h*(k+1))<=t:
            d+=(Find_G(k*h+gamma*h, k1)-Find_G(k*h, k1))/E1
            d+=(Find_G(k*h+h, k1)-Find_G(k*h+gamma*h, k1))/E2
            if (h*(k+1))==t:
                break
        else:
            if h*(k)+gamma*h>=t:
                d+=(Find_G(t, k1)-Find_G(k*h, k1))/E1
            else:
                d+=(Find_G(k*h+gamma*h, k1)-Find_G(k*h, k1))/E1
                d+=(Find_G(t, k1)-Find_G(k*h+gamma*h, k1))/E2
            break
    return d
def displacement_exact_fast(t, k1):
    d=0
    d1=0
    d2=0
    if k1 == 0:
        return 0
    for k in range (0,int(t//h)):
        d+=(Find_G(k*h+gamma*h, k1)-Find_G(k*h, k1))/E1
        d+=(Find_G(k*h+h, k1)-Find_G(k*h+gamma*h, k1))/E2
    k=int(t//h)
    if h*(t//h)+gamma*h>=t:
        d+=(Find_G(t, k1)-Find_G(k*h, k1))/E1
    else:
        d+=(Find_G(k*h+gamma*h, k1)-Find_G(k*h, k1))/E1
        d+=(Find_G(t, k1)-Find_G(k*h+gamma*h, k1))/E2
    return d
def displacement_deg_1(t):
    return find_w0(t)+h*(find_w1(t)+N1((t%h)/h)*find_w0_der(t))
def displacement_deg_2(t):
    return displacement_deg_1(t)+h**2*(find_w2(t)+N1((t%h)/h)*find_w1_der(t)+N2((t%h)/h)*find_w0_der2(t))

def find_w0(t):
    return b*(b*(math.log(4*b*b+(L-2*t)**2)-math.log(4*b*b+L*L))-(L-2*t)*math.atan((L-2*t)/(2*b))+L*math.atan(L/(2*b)))/h0
def find_w0_der(t):
    return (find_w0(t+L/points)-find_w0(t-L/points))/(2*L/points)
def find_w0_der2(t):
    return (find_w0_der(t+L/points)-find_w0_der(t-L/points))/(2*L/points)
def find_w1(t):
    return -N1(0)*find_w0_der(0)+t*2/L*N1(0)*find_w0_der(0)
def find_w1_der(t):
    return (find_w1(t+L/points)-find_w1(t-L/points))/(2*L/points)
def find_w2(t):
    return -N1(0)*find_w1_der(0)-N2(0)*find_w0_der2(0)

def N1(t):
    if t <= gamma:
        return h0*t/E1-t+c2
    else:
        return h0*(gamma/E1+(t-gamma)/E2)-t+c2
def N2(t):
#   N1 = Ax+B | x < gamma,
#   N1 = A*gamma+B+C*(x-gamma) | x >= gamma
    A = h0/E1-1; B = c2; C = h0/E2-1; D = A*gamma+B-gamma*C
    c3 = A/3*gamma**3+B/2*gamma**2+C/3+D/2-C/3*gamma**3-D/2*gamma**2
    if t <= gamma:
        return A/2*t**2+B*t-c3
    else:
        return A/2*gamma**2+B*gamma-c3+C/2*t**2+D*t-c3-C/2*gamma**2-D*gamma

def Find_stress(t):
    pass
if __name__ == '__main__':
    pool = Pool()
    start_time = time.time()
    Find_k1()
    print (f"k1 = {k1}")
    fig = plt.subplots()
    plt.subplot(2, 1, 1)
    x_local = [0.5-1/2/scale+i/points/scale for i in range(0,points+1)]
    x_global = [i/points for i in range(0,points+1)]
    #y = [displacement_exact(x_local[i]) for i in range (0,points+1)]
    #d1 = [displacement_deg_1(x_local[i]) for i in range (0,points+1)]
    d2 = [displacement_deg_2(x_local[i]) for i in range (0,points+1)]
    print("th start")
    t_th_start=time.time()
    d_e_th = list(pool.starmap(displacement_exact_fast, zip(x_local, [k1 for i in range(len(x_local))])))
    pool.close()
    pool.join()
    pool = Pool()
    d1_th = list(pool.map(displacement_deg_1, x_local))
    pool.close()
    pool.join()
    pool = Pool()
    d2_th = list(pool.map(displacement_deg_2, x_local))
    pool.close()
    pool.join()
    print("th time  --- %s seconds ---" % (time.time() - t_th_start))
    #plt.plot(x_local, y, "y")
    plt.plot(x_local, d_e_th, "y")
    plt.plot(x_local, d1_th, "b--")
    plt.plot(x_local, d2_th, "r--")
    plt.subplot(2, 1, 2)
    Force_arr = list(map(Find_force, x_global))
    plt.plot(x_global, Force_arr)
    #plt.plot(x, w0)
    #plt.plot(x_local, d1, "b--")
    #plt.plot(x_local, d2, "r--")
    #N1_arr = [N1(i/points) for i in range (0,points+1)]
    #N2_arr = [N2(i/points) for i in range (0,points+1)]
    #print(f"c2 = {c2}")
    #plt.plot(x_local, N2_arr)
    #plt.plot(x, N1_arr)
    print("Global   --- %s seconds ---" % (time.time() - start_time))
    plt.show()
