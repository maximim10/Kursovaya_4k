import math
import matplotlib.pyplot as plt
import time
L=1; n=200; h=L/n; b=h*1; gamma=1/2; E1=2; E2=E1/10; points=100000; scale=n
h0 = 1/(gamma/(E1)+(1-gamma)/(E2)); c2 = 1/2-h0*(gamma**2/E1/2+(1-gamma)*gamma/E1+(1-gamma)**2/E2/2)
start_time = time.time()
a = h0/E1-1; b = c2; c = h0/E2-1; d = a*gamma+b-gamma*c
c3 = a/3*gamma**3+b/2*gamma**2+c/3+d/2-c/3*gamma**3-d/2*gamma**2
k1=0
def Find_G_part(t):
    G=b**2*math.log(4*b**2+(L-2*t)**2)-b*(L-2*t)*math.atan((L-2*t)/(2*b))
    return G
def Find_G(t):
    G=Find_G_part(t)+k1*t
    return G
def Find_k1():
    global k1
    summ=0
    for k in range (0,n):
        summ+=(Find_G_part(k*h+gamma*h)-Find_G_part(k*h))/E1
        summ+=(Find_G_part(k*h+h)-Find_G_part(k*h+gamma*h))/E2
    k1=-summ/Sum_near_k1()
    print (f"k1 = {k1}")
    #return k1
def Sum_near_k1():
    summ=0
    for k in range (0,n):
        summ+=gamma*h/E1
        summ+=h*(1-gamma)/E2
    return summ


def displacement_exact(t):
    d=0
    for k in range (0,n):
        if (h*(k+1))<=t:
            d+=(Find_G(k*h+gamma*h)-Find_G(k*h))/E1
            d+=(Find_G(k*h+h)-Find_G(k*h+gamma*h))/E2
            if (h*(k+1))==t:
                break
        else:
            if h*(k)+gamma*h>=t:
                d+=(Find_G(t)-Find_G(k*h))/E1
            else:
                d+=(Find_G(k*h+gamma*h)-Find_G(k*h))/E1
                d+=(Find_G(t)-Find_G(k*h+gamma*h))/E2
            break
    return d
def displacement_deg_1(t):
    return find_w0(t)+h*(find_w1(t)+N1((t%h)/h)*find_w0_der(t))
def displacement_deg_2(t):
    return displacement_deg_1(t)+h**2*(find_w2(t)+N1((t%h)/h)*find_w1_der(t)+N2((t%h)/h)*find_w0_der2(t))

def find_w0(t):
    return b*(b*(math.log(4*b*b+(L-2*t)**2)-math.log(4*b*b+L*L))-(L-2*t)*math.atan((L-2*t)/(2*b))+L*math.atan(L/(2*b)))/h0
def find_w0_der(t):
    #if not t == 0 and not t == L:
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
#   N1 = ax+q | x < gamma,
#   N1 = a*gamma+q+c*(x-gamma) | x >= gamma
    
    if t <= gamma:
        return a/2*t**2+b*t-c3
    else:
        return a/2*gamma**2+q*gamma-c3+c/2*t**2+d*t-c3-c/2*gamma**2-d*gamma
print("1--- %s seconds ---" % (time.time() - start_time))        
Find_k1(); print("2--- %s seconds ---" % (time.time() - start_time))
fig = plt.subplots()
x_local = [0.5-1/2/scale+i/points/scale for i in range(0,points+1)]
y = [displacement_exact(x_local[i]) for i in range (0,points+1)]; print("3--- %s seconds ---" % (time.time() - start_time))
d1 = [displacement_deg_1(x_local[i]) for i in range (0,points+1)]; print("4--- %s seconds ---" % (time.time() - start_time))
d2 = [displacement_deg_2(x_local[i]) for i in range (0,points+1)]; print("5--- %s seconds ---" % (time.time() - start_time))
#w0 = [find_w0(x_local[i]) for i in range (0,points+1)]
#w1 = [find_w1(x_local[i]) for i in range (0,points+1)]

plt.plot(x_local, y, "y") 
#plt.plot(x, w0)
plt.plot(x_local, d1, "b--")
plt.plot(x_local, d2, "r--")
#N1_arr = [N1(i/points) for i in range (0,points+1)]
#N2_arr = [N2(i/points) for i in range (0,points+1)]
#print(f"c2 = {c2}")
#plt.plot(x, N2_arr)
#plt.plot(x, N1_arr)
print("6--- %s seconds ---" % (time.time() - start_time))
plt.show()
