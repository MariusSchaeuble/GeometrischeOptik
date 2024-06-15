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
sigma_b = 1/100
g = 10.5/100
sigma_g = 1/100

f1 = b*g/(b+g)
sigma_f1 = gauss("b*g/(b+g)")

b = 10.3/100
g = (30-16.6)/100
sigma_b = 1/100
sigma_g = 1/100

f2 = b*g/(b+g)
sigma_f2 = gauss("b*g/(b+g)")

b = (74-30)/100
g = (30-23.5)/100
sigma_b = 1/100
sigma_g = 1/100

f3 = b*g/(b+g)
sigma_f3 = gauss("b*g/(b+g)")

b = (38.7-30)/100
g = (30-11.5)/100
sigma_b = 1/100
sigma_g = 1/100

f4 = b*g/(b+g)
sigma_f4 = gauss("b*g/(b+g)")

#system aus 60, -150 mm
b = (47.1-30)/100
g = (30-11.5)/100
sigma_b = 1/100
sigma_g = 1/100
f1s = b*g/(b+g)
sigma_f1s = gauss("b*g/(b+g)")

b=(53.9-30)/100
g=(30-17.8)/100
sigma_b = 1/100
sigma_g = 1/100

f2s = b*g/(b+g)
sigma_f2s = gauss("b*g/(b+g)")

b=(42.8-30)/100
g=(30-7.9)/100
sigma_b = 1/100
sigma_g = 1/100

f3s = b*g/(b+g)
sigma_f3s = gauss("b*g/(b+g)")

b=(67.1-30)/100
g=(30-19.6)/100
sigma_b = 1/100
sigma_g = 1/100

f4s = b*g/(b+g)
sigma_f4s = gauss("b*g/(b+g)")

#3.1.2 teil: Bessel
#abstande g+b=s
#System aus 60mm, -150mm
s=59/100
e=38/100
sigma_e = 1/100
sigma_s = 1/100

f1sb = (s**2 - e**2)/(4*s)
sigma_f1sb = gauss("(s**2 - e**2)/(4*s)")

s=52.1/100
e=30.7/100

f2sb = (s**2 - e**2)/(4*s)
sigma_f2sb = gauss("(s**2 - e**2)/(4*s)")

s=39/100
e=14.5/100

f3sb = (s**2 - e**2)/(4*s)
sigma_f3sb = gauss("(s**2 - e**2)/(4*s)")

s=68.9/100
e=(76.6-27.4)/100

f4sb = (s**2 - e**2)/(4*s)
sigma_f4sb = gauss("(s**2 - e**2)/(4*s)")

#sammellinse, 100mm
s=(68.9)/100
e=(73.1-27.4)/100

f1b = (s**2 - e**2)/(4*s)
sigma_f1b = gauss("(s**2 - e**2)/(4*s)")

s=(60)/100
e=(63.5-28)/100

f2b = (s**2 - e**2)/(4*s)
sigma_f2b = gauss("(s**2 - e**2)/(4*s)")

s=(51.7)/100
e=(54.1-28.7)/100

f3b = (s**2 - e**2)/(4*s)
sigma_f3b = gauss("(s**2 - e**2)/(4*s)")

s=(41.7)/100
e=(42-31.3)/100

f4b = (s**2 - e**2)/(4*s)
sigma_f4b = gauss("(s**2 - e**2)/(4*s)")

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
s0 = 250
f2 = 20
sigma_f2 = 2
sigma_s0 = 25
f1 = 40
sigma_f1 = 4

V1 = 2.2*s0/f2
sigma_V1 = gauss("2.2*s0/f2")

l = 228
sigma_l = 10

V2 = s0*(l-f2)/(f1*f2)
sigma_V2 = gauss("s0*(l-f2)/(f1*f2)")

V3 = s0*l/(f1*f2)
sigma_V3 = gauss("s0*(l)/(f1*f2)")

