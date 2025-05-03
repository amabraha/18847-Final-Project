import math
N = 20
points = []

for theta in [math.pi*i/(2*N) for i in range(N+1)]:
    points += [[0+5*math.cos(theta),9+5*math.sin(theta)]]
    
for theta in reversed([math.pi*i/(2*N) for i in range(N+1)]):
    points += [[0+6*math.cos(theta),9+6*math.sin(theta)]]

M = len(points)
print(str(M),"2 0 0")
for i in range(1,M+1):
    print(str(i),points[i-1][0],points[i-1][1])

print(str(M),"0")
print("1",str(M),"1")
for i in range(2,M+1):
    print(str(i),str(i-1),str(i))
print("0\n")
