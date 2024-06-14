import matplotlib
import numpy
import numpy as np
import sympy as sym
from Helpers import identifier, isCharacter
import math
from numpy import matrix, array, mean, std, max, linspace, ones, sin, cos, tan, arctan, pi, sqrt, exp, arcsin, arccos, arctan2, sinh, cosh
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show, xlabel, ylabel, legend, title, savefig, errorbar, grid
import scipy.optimize as opt
from GPII import *
from math import sqrt
pi = math.pi


matplotlib.rc('xtick', labelsize=20)
matplotlib.rc('ytick', labelsize=20)







def gauss(term):
    ids = identifier(term)
    symbols = []
    for str1 in ids:
        symbols.append(sym.sympify(str1))
    termSymbol = sym.sympify(term)
    values = []
    for ident in ids:
        exec("values.append(" + ident + ")")

    derivatives = []
    i = 0
    while i < len(symbols):
        r = sym.diff(termSymbol, symbols[i])
        j = 0
        while j < len(symbols):
            # exec('r.evalf(subs={symbols[j]: ' + values[j] + '})')
            r = r.evalf(subs={symbols[j]: values[j]})
            j += 1
        derivatives.append(r.evalf())
        i += 1
    i = 0
    while i < len(derivatives):
        exec("derivatives[i] *= sigma_" + ids[i])
        i = i + 1
    res = 0
    for z in derivatives:
        res += z ** 2
    return math.sqrt(res)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#3.1.1 teil: abbildungsgleichung

#sammellinse, 60mm
b = 15.4/100
g = 10.5/100

b = 10.3/100
g = (30-16.6)/100

b = (74-30)/100
g = (30-23.5)/100

b = (38.7-30)/100
g = (30-11.5)/100

#system aus 60, -150 mm
b = (47.1-30)/100
g = (30-11.5)/100

b=(53.9-30)/100
g=(30-17.8)/100

b=(42.8-30)/100
g=(30-7.9)/100

b=(67.1-30)/100
g=(30-19.6)/100

#3.1.2 teil: Bessel
#abstande g+b=s
#System aus 60mm, -150mm
s=59/100
e=38/100

s=52.1/100
e=30.7/100

s=39/100
e=14.5/100

s=68.9/100
e=(76.6-27.4)/100

#sammellinse, 100mm
s=(68.9)/100
e=(73.1-27.4)/100

s=(60)/100
e=(63.5-28)/100

s=(51.7)/100
e=(54.1-28.7)/100

s=(41.7)/100
e=(42-31.3)/100

#3.1.3 teil: Autokollimation (systematischer fehler: lochblende um 1cm verschoben zur messapperatur, daher(x+1))
f_100 = (23.9-(12.9+1))/100
f_60 = (19.7-(12.8+1))/100
f_sys = (23.3-(12.8+1))/100

#3.2.1
#system: Mattscheibe -> blende 2mm -> Sammellinse, 40mm: strahlaufbreiterung nicht zu vermeiden
#3.2.2
#system: Mattscheibe -> blende 2mm-> Sammellinse, 40mm -> sammellinse 60mm -> lochblende 0.6mm -> 80mm: strahlverbreiterung aufgehoben, punkt weiterhin unscharf, aber bei jedem abstand

#3.3
#f_objektiv= 100mm, f_okular=1.5mm

#3.4
#zwischenbild mit skala vergleich: 2.2 fache vergrößerung
#objektiv: 40mm
#okular: 20mm
#Abstand der Linsen: 22.8cm