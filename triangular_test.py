import numpy as np
import scipy as sp

Lx = 30 ; Ly = 30
n = Lx * Ly
t = 1.0

H = np.zeros ((n,n))

for iy in range (Ly):
  for ix in range (Lx-1):
    i = iy*Lx + ix
    jx = ix + 1
    if jx >= Lx:
      jx = 0
    j = iy*Lx + jx
    H[i,j] = -t
    H[j,i] = -t
    
for iy in range (Ly-1):
  for ix in range (Lx):
    i = iy*Lx + ix
    jy = iy + 1
    if jy >= Ly:
      jy = 0
    j = jy*Lx + ix
    H[i,j] = -t
    H[j,i] = -t
    
for iy in range (Ly-1):
  for ix in range (Lx-1):
    i = iy*Lx + ix
    jx = ix + 1
    jy = iy + 1
    if jx >= Lx:
      jx = 0
    if jy >= Ly:
      jy = 0
    j = jy*Lx + jx
    H[i,j] = -t
    H[j,i] = -t
    
W = sp.linalg.eigvalsh (H)

# for j in range (n):
#   for i in range (n):
#     if (H[i,j] != 0.0):
#       print(i+1, '    ',j+1, '    ', H[i,j])

with open('eigenvals_test.txt', "w") as file:
  for i in range (len(W)):
    file.write('%6d %14.6f \n' %(i, W[i]))