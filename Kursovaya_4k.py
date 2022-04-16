import math
import matplotlib.pyplot as plt
import time
L=1; n=20; h=L/n; b=5*h; gamma=1/2; E1=2; E2=E1/2; points=10000; scale=n/4
h0=1/(gamma/(E1)+(1-gamma)/(E2))
start_time = time.time()
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
def displacement(t):
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
def find_w0(t):
    return b*(b*(math.log(4*b*b+(L-2*t)**2)-math.log(4*b*b+L*L))-(L-2*t)*math.atan((L-2*t)/(2*b))+L*math.atan(L/(2*b)))/h0
def find_w0_der(t):
    #if not t == 0 and not t == L:
    return (find_w0(t+L/points)-find_w0(t-L/points))/(2*L/points)
def find_w0_der2(t):
    return (find_w0_der(t+L/points)-find_w0_der(t-L/points))/(2*L/points)
def find_w1(t):
    return -N1(0)*find_w0_der(0)+t*2/L*N1(0)*find_w0_der(0)
def N1(t):
    if t <= gamma:
        return h0*t/E1-t+
Find_k1()
fig = plt.subplots()
x = [i/1000 for i in range(0,1001)]
y = [displacement(i/1000) for i in range (0,1001)]
w0 = [find_w0(i/1000) for i in range (0,1001)]
w1 = [find_w1(i/1000) for i in range (0,1001)]
plt.plot(x, y)
plt.plot(x, w0)
plt.plot(x, w1)
plt.show()
print("--- %s seconds ---" % (time.time() - start_time))
