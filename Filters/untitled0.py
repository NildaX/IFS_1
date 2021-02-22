# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 18:21:07 2020

@author: hp
"""
import numpy as np
import os
#import matplotlib
import matplotlib.pyplot as plt

nombre3 = 'i-SDSS-III.txt'
print(nombre3)
lineas = []
if os.path.isfile(nombre3):
    print('El archivo existe.');
    archivo = open(nombre3, 'r')
    for linea in archivo.readlines():
        if linea.startswith("#")==False and linea.isalnum()==False and linea!= '':
            linea = linea.rstrip('\r\n')
            lineas.append(linea)
            
   #     print(lineas)
    lineas = list(filter(bool,lineas))
    arr_w = np.zeros(len(lineas))
    arr_f =  np.zeros(len(lineas))
    j= 0
    for i in lineas:
        primeralinea=i.split()
        arr_w[j] = primeralinea[0]
        arr_f[j] = primeralinea[1]
        j = j+1
            
    print(arr_w)
    print(arr_f)
    plt.plot(arr_w,arr_f,'--')
    plt.show()
    archivo.close()
else:
    print("el archivo no existe")