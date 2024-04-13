import numpy as np
import matplotlib.pyplot as plt

n = 30 #number of items
la = 10 #average number of arrivals in unit time interval
m_station = 4

td = np.random.exponential(1 / la, n) #Exp(la)-distributed arrival time differences
#print(td)
#service times, uniformly distributed between d-delta ... d+delta
d=0.25
delta=d/5
m=d-delta
M=d+delta
ts=np.random.uniform(m,M,n)

ta=np.zeros(n) #arrival times
tsc=np.zeros(n) #service completed times

#aa=np.zeros(n) #arrival times
#bb=np.zeros(n) #service completed times
service_starting_time = np.zeros(n)

tsc[0]=ta[0]+ts[0]
"""
bb[0] = aa[0] + ts[0]

for k in range(n-1):
    aa[k+1]=aa[k]+td[k]
    bb[k+1]=bb[k+1]+ts[k+1]
"""


# station
station = np.zeros(m_station)
color = np.zeros(n) # at the index of k (1 -> n) what color it should be
cnt = 0 # count at what color is it

# waiting
waiting_start = np.zeros(n)
waiting_end = np.zeros(n)

color[0] = 1
color[1] = 2
color[2] = 3
color[3] = 4
total_time = 0
for k in range(n):
  total_time += td[k] + waiting_end[k] - waiting_start[k]
total_time = round(total_time, 2)

for k in range(0,3):
  ta[k + 1] = ta[k] + td[k]
  tsc[k + 1] = ta[k + 1] + ts[k + 1]
  station[k] = tsc[k]

station[3] = tsc[3]
# station
# 0 1 2 3


for k in range(3, n-1):
    ta[k + 1] = ta[k] + td[k]

    nk = np.argmin(station)

    if(station[nk] > ta[k + 1]):
      waiting_start[k + 1] = ta[k + 1]
      waiting_end[k + 1] = station[nk]
      ta[k + 1] = station[nk]

    color[k + 1] = nk + 1
    tsc[k + 1] = ta[k + 1] + ts[k+1] + waiting_end[k + 1] - waiting_start[k + 1]
    station[nk] = tsc[k + 1]

# number of items in service
tasc=np.hstack((ta,tsc)) #arrival and service completed times
asc=np.hstack((np.zeros(n),np.ones(n))) #0 = arrives, 1 = service completed
tascs=np.sort(tasc) #smallest->largest
ind=np.argsort(tasc) #corresponding indices
asct=asc[ind] #arrivals and service completions in time order
nis=np.zeros(2*n-1) #number of items in service between asct time-intervals
nisp=0 #previous number of items in service

sum = 0
for k in range(2*n-1):
    if asct[k]==0:
        nis[k]=nisp+1 #arrival
    else:
        nis[k]=nisp-1 #service completed
    nisp=nis[k]

for k in range(2*n-1):
  sum += nis[k]
sum = round(sum / 60, 2)

# ta la diem bat dau, tsc diem ket thuc, td random thoi gian bat dau, ts random thoi gian ket thuc
plt.figure(figsize=(15,10))
plt.subplot(311)

plt.plot([ta[0],tsc[0]],[0,0], 'r', label = 'station 1' ,lw=3)

plt.plot([ta[1],tsc[1]],[1,1],'g', label = 'station 2' ,lw=3)
plt.plot([ta[2],tsc[2]],[2,2],'b', label = 'station 3' ,lw=3)
plt.plot([ta[3],tsc[3]],[3,3],'purple', label = 'station 4'  ,lw=3)

ok = False
for k in range(4, n):
    if(color[k] == 1):
      plt.plot([ta[k],tsc[k]],[k,k],'r' ,lw=3)
    elif(color[k] == 2):
      plt.plot([ta[k],tsc[k]],[k,k],'g' ,lw=3)
    elif(color[k] == 3):
      plt.plot([ta[k],tsc[k]],[k,k],'b' ,lw=3)
    elif(color[k] == 4):
      plt.plot([ta[k],tsc[k]],[k,k],'purple' ,lw=3)

    if(waiting_start[k] != 0 and ok == False):
        ok = True
        plt.plot([waiting_start[k],waiting_end[k]],[k,k],'black', label = 'waiting' ,lw=3)
    if(ok == True):
      plt.plot([waiting_start[k],waiting_end[k]],[k,k],'black' ,lw=3)
# k is the hang ngang

plt.plot(ta,-np.ones(n),'r.',label='arrives')
plt.plot(tsc,-np.ones(n),'b.',label='service completed')


plt.xlim(0,tascs[-1])
plt.grid()
plt.legend(fontsize=14)
plt.title('n = '+ str(n) + ', ' +'$\lambda = $'+str(la)+r', $d$ = '+str(d)
          +', $\delta$ = '
          +str(delta) + ', m = ' + str(m_station) + ': total time ' + str(total_time) ,fontsize=16)
plt.yticks(np.arange(0,n+1,5))
plt.ylabel('item number',fontsize=16)


plt.subplot(312)

plt.bar(tascs[:-1],nis,zorder=2,facecolor='b',
        align='edge', #tolppien vasen reuna tsps[:-1]:ssä
        width=tascs[1:]-tascs[:-1]) #tolppien leveydet


plt.plot(ta,-0.1+np.zeros(n),'r.',label='start')
plt.plot(tsc,-0.1+np.zeros(n),'b.',label='end')
plt.grid()
plt.title('average number of items being serviced: ' + str(sum))
plt.xlim(0,tascs[-1])
plt.yticks(np.arange(0,np.max(nis)+1,1))
plt.ylabel('items in service',fontsize=16)


ss=np.hstack((waiting_start,waiting_end)) #arrival and service completed times
uu=np.hstack((np.zeros(n),np.ones(n))) #0 = arrives, 1 = service completed
sss=np.sort(ss) #smallest->largest
dd=np.argsort(ss) #corresponding ddices
uuu=uu[dd] #arrivals and service completions in time order
zz=np.zeros(2*n-1) #number of items in service between asct time-intervals
z=0 #previous number of items in service

sum = 0
for k in range(2*n-1):
    if uuu[k]==0:
        zz[k]=z+1 #arrival
    else:
        zz[k]=z-1 #service completed
    z=zz[k]

for k in range(2*n-1):
  if(zz[k] < 0):
    zz[k] = 0
  sum += zz[k]

plt.subplot(313)

sum = round(sum / 60, 2)

plt.bar(sss[:-1],zz,zorder=2,facecolor='black',
        align='edge', #tolppien vasen reuna tsps[:-1]:ssä
        width=sss[1:]-sss[:-1]) #tolppien leveydet

plt.xlim(0,tascs[-1])
plt.plot(ta,-0.1+np.zeros(n),'r.',label='start')
plt.plot(tsc,-0.1+np.zeros(n),'b.',label='end')
plt.grid()
plt.title('average number of items waiting: ' + str(sum))
plt.yticks(np.arange(0,np.max(nis)+1,1))
plt.ylabel('items waiting',fontsize=16)
plt.xlabel('time t',fontsize=16)
plt.show()

