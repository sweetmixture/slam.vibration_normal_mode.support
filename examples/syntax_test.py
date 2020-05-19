#!/bin/python

import numpy as np

arr = np.matrix([[1,-1],[1,1]])

arr2 = np.transpose(arr)

print(arr)
print(arr2)



mat = np.matrix( [[1,2,2],[-1,0,2],[0,0,0]] )
print(mat)

q, r = np.linalg.qr( mat )

print(q)
print(r)
