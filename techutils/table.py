# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 19:20:35 2018

@author: Avedis Neehus
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from statsmodels.stats.weightstats import DescrStatsW as data
from uncertainties import ufloat as u, umath
from pint import  UnitRegistry  
import sympy as sym
import inspect

ur = UnitRegistry()
Q_ = ur.Quantity
Q_.pm = Q_.plus_minus
Q_.mag= Q_.magnitude 

class attriList(list):
    
    def __getattr__(self,attr):
        
            return attriList(item.__getattribute__(attr) for item in self)
        
def listify(f):
        
    def wrapping(*args,**kwargs):        
        return attriList(map(lambda *arg : f(*arg,**kwargs), *args))
    
    return wrapping

def texify( *varnames,fname = None, **cnames):    
    def deco(f):        
        def function(*args):
            
            entity  = fname if fname else f.__name__
            
            argspec =  inspect.getfullargspec(f)            
            def_args = argspec.args[:-len(argspec.defaults)] if argspec.defaults else argspec.args
            def_kwargs = argspec.args[-len(argspec.defaults):] if argspec.defaults else {}
            
            
            symargs   = tuple(sym.Symbol(varname)  for varname in varnames)  if varnames else tuple(sym.symbols(def_args))
            symkwargs = tuple(sym.Symbol(cname) for cname in cnames.values())if cnames   else tuple(sym.symbols(def_kwargs))
            symall = symargs+symkwargs
            symexpr = f(*symall)
            
            allarg = args+f.__defaults__ if f.__defaults__ else args
            
            
            numargs = tuple(sym.Symbol('{:Lx}'.format(arg.value)) if  hasattr(arg, 's') else sym.Symbol('{:Lx}'.format(arg)) for arg in allarg  )
            numexpr = f(*numargs)
            
            print('\\begin{align*}')
            print(entity + ' &= '+ sym.latex(symexpr) + r'\\')
            print('&=' + sym.latex(numexpr, mul_symbol = 'dot') + '=', end = '')
            try:
                result = f(*args)
                print( '{:Lx.02f}'.format(result.value))
            except: pass
            print('\\end{align*}')
            
        return function
    return deco

def techutil(*varnames,fname = None, **cnames):
    
    def deco(f):
                
        def function(*args):            
            return f(*args)    
        
        function.gauß = errify(*varnames,fname = fname, **cnames)(f)
        function.tex = texify(*varnames,fname= fname , **cnames)(f)   
        
        return function    
    return deco


def errify( *varnames, fname = None, **cnames):
    
    def deco(f):
        
        def function(*args):
            
            entity  = fname if fname else f.__name__
            
            argspec =  inspect.getfullargspec(f)            
            def_args = argspec.args[:-len(argspec.defaults)] if argspec.defaults else argspec.args
            def_kwargs = argspec.args[-len(argspec.defaults):] if argspec.defaults else {}
            
            
            symargs   = tuple(sym.Symbol(varname)  for varname in varnames)  if varnames else tuple(sym.symbols(def_args))
            symkwargs = tuple(sym.Symbol(cname) for cname in cnames.values())if cnames   else tuple(sym.symbols(def_kwargs))
            symall = symargs+symkwargs
            
            allarg = args+f.__defaults__ if f.__defaults__ else args
            
            
            symprecise = [symbol for symbol,arg in zip(symall,allarg) if hasattr(arg, 's')]
            expr = f(*symall)#, **symkwargs)
            
            derivatives = [expr.diff(symbol) for symbol in symprecise]
            
            
            inprecise = [sym.Symbol('{:Lx}'.format(arg.error)) for arg in allarg if hasattr(arg, 's')]
            numargs = tuple(sym.Symbol('{:Lx}'.format(arg.value)) if  hasattr(arg, 's') else sym.Symbol('{:Lx}'.format(arg)) for arg in allarg  )
            numexpr = f(*numargs)
            numderivs =  [numexpr.diff(symbol) for symbol,val in zip(numargs,allarg) if hasattr(val, 's')]
            
            
            print('Der Fehler von ' + entity+ ' ist:')
            print('\\begin{align*}')
            
            print('\\Delta ' + entity + '&= \\sqrt{', end =  '')
            
            for symbol in symprecise:
                print('\\left(' + '\\frac{\\partial ' + entity + '}{\\partial ' + sym.latex(symbol) + '}'  + '\\Delta ' + sym.latex(symbol) +  ' \\right)^{2}' + '+', end = ' ' ) 
                            
            print(r' } \\ &= ')
            
            print(r'\sqrt{', end = '')
            for val,derivativ in zip(symprecise,derivatives):
                print('\\left(' + sym.latex(derivativ) + '\\right)^{2}' + ' \\cdot ' + '\\left( \\Delta ' + sym.latex(val) + '\\right)^{2} + ', end = ' ' )
            print(r'} \\')    
            
            print(r'&= \sqrt{', end = '')
            for val,numderiv in zip(inprecise,numderivs):
                print('\\left(' + sym.latex(numderiv,mul_symbol='dot') + '\\right)^{2}' + ' \\cdot ' + '\\left(' + sym.latex(val) + '\\right)^{2} +', end = ' ' )
            print(r'} \\')
            print('\\Delta ' + entity )
            
             
            
            try:        
                value = f(*args)
                print(' &= ' , '{:Lx.03f}'.format(value.s * value.units))  
            except:
                pass
            print(r'\end{align*}')
            
        return function 
    return   deco          
        
        
        

def add(a, b):    
    res = []
    
    if type(a[0]) == list:         
        for v in a: 
            
            res += list(map(lambda x: x+b, v))
    else: 
        res = list(map(lambda x: x+b, a))
      
    return np.array(res)



def runde(a, b  = 2):
    
    res = []
    
    if type(a[0]) == list:    
        for v in a:        
            res+= list(list(map(lambda x: round(x,b), v)))            
    else:
        res =  list(map(lambda x: round(x,b), a))
        
    return res
    
#class stats(data):
    
 #   @property
#    def sigma(self, dig = 2):
 #       return runde(self.std/np.sqrt(self.sum_weights), dig)
  



def plotit(data,xaxis, xlabel = None, ylabel = None, maxy= None, labels = None):
    
    if labels:
        for seq,lab in zip(data, labels):
            plt.bar(xaxis, seq, alpha=0.6, label =lab, edgecolor = 'b')
    else:
        for seq in data:
            plt.bar(xaxis, seq, alpha=0.6, edgecolor = 'b')
            
    if xlabel: 
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if maxy:
        plt.yticks(np.arange(maxy))
    plt.show()    
    
    #print(latex_figure(savefig()))

def tabelize(data, collab = None, rowlab = None, hline = 0):
    
    form = '{' + ''.join('c' for a in data[0])
    
    if not(rowlab is None):
        form += 'c}'
    else:
        form = form.join(' }')
                    
      
    print(r'\begin{tabular} ' + form )
    
    if not(collab is None):
        clabels = r''.join(' & ' + a for a in collab)
        print(clabels+ r'\\')
    
    if not(rowlab is None): 
        
        for vals,rlabel in zip(data,rowlab):
            row = str(rlabel) + ''.join(' & ' + str(val) for val in vals)
                        
            print(row+r'\\')
            if hline:
                print(r'\hline')
    else: 
        
        for vals in data:
            row = r''.join(str(val) + ' & ' for val in vals)  
            row = row[:-2]
            print(row + r'\\') 
            if hline:
                print(r'\hline')
      
    print(r'\end{tabular}')

@listify
def mean(li, dig = 1):
    try:
        val = round(np.mean(li), dig)
        
    except TypeError:
        
        val = np.mean(li)
        
    if dig:
        return val    
    else:
        return int(val)
    

def std(li, dig = 1):
    
    val = round(np.std(li), dig)

    if dig:
        return val
    
    else:
        return int(val)
    
def sigma(li, n = None, dig = 1 ):
    
    if n:
        val = np.std(li, ddof = 1)/np.sqrt(n)
    else:
        val = np.std(li, ddof = 1)/np.sqrt(len(li))
        
    if dig:
        
        try:
            
            return round(val,dig)
        
        except TypeError:
            return val
    
    else:
        return int(val)
    
    class mpi(float):
        
       def  __mul__(self, other):
#           if not(type(other :
               return self * other
          # else: 

def mprint(q):
    print( '{:Lx}'.format(q))

def ChiSquare(f,xdata,ydata,popt):
    res = np.sum((ydata-f(xdata,*popt))**2)
    return res/(len(ydata)-len(popt))

def Rsquare(f,xdata,ydata,popt):
    res = np.sum((ydata-f(xdata,*popt))**2)
    tot = np.sum((ydata-np.mean(ydata))**2)
    return 1-res/tot
def SE(f,xdata,ydata,popt):
    res = np.sum((ydata-f(xdata,*popt))**2)
    return (res/(len(ydata)))**0.5
def parerror(cov):
    return np.sqrt(np.diag(cov))

T1 = np.array([[878,650,550,575,653,765,907],
               [888,644,550,569,647,769,897],
               [890,644,547,578,659,759,900]])*10**-3
    
T2 = np.array([[958,831,700,594,534,534,597,694],
               [971,831,703,597,538,534,594,697],
               [972,842,694,597,541,534,590,690]])*10**-3
    
T1 = (T1[0]+T1[1]+T1[2] )/3
T1 = T1[2:]
l1 = np.arange(0,5)
l2 = np.arange(-7,1)
mz = Q_('172.2g').pm(0.1)
mk = Q_('182.2g').pm(0.1)

def f(x,a,b):
    return a + b*x

popt, pcov = curve_fit(f, l1**2, T1**2,absolute_sigma = 1, sigma = [2*t*0.006 for t in T1] )
print(Rsquare(f, l1**2, T1**2, popt),ChiSquare(f, l1**2, T1**2, popt) )
#plt.errorbar(l1**2,T1**2, yerr = [2*t*0.006 for t in T1], fmt = 'g3', capsize = 3)   
l1 = np.linspace(0,17) 
#plt.plot(l1, f(l1,*popt))
plt.xlabel('$a$ in cm ')
plt.ylabel('$T$ in s')
    
b =  Q_(str(popt[1])+'s**2/cm**2').pm(pcov[1][1]**0.5)
a =  Q_(str(popt[0])+'s**2').pm(pcov[0][0]**0.5)
@techutil(fname = 'D^{*}',b = 'B',mz = 'm_z')
def D(b = b, mz = mz):
    return 4*np.pi**2*mz*1/b

@techutil(fname = 'I_z',r ='R',mz ='m_z')
def iz(r=Q_('5.19/2cm').pm(0.05), mz = mz):
    return 1/2 * mz * r**2

@techutil(fname ='I_t', a = 'A', d = 'D', iz = 'I_z')
def It(a  =  a,d = D().to('m**2*g*s**-2'),iz = iz().to('g*m**2')):
    return (1/(4*np.pi**2))*a*d-iz

@techutil(fname = 'I_p',T = 'T_{min}', D ='D^{*}', it = 'I_t' )
def minna(T = Q_('0.529s').pm(0.004) ,D = D().to('m**2*g*s**-2'),it = It()):
    return T**2*(D/(4*np.pi**2))-it


def f2(l,ik,):
    l = 0.01*(l+2.5057)*ur.meter
    ik = ik*It().units
    return (2*np.pi*np.sqrt(1/D().value.to_base_units() *(ik.to_base_units()+It().value.to_base_units()+mk.value.to_base_units()*l**2))).magnitude 

def f3(l,a,b,c,d):
   return a*(l-b)**2+ c*(l-b)**4 +d

def f4(l,ik,o, m):
    l = 0.01*(l-o)*ur.meter
    ik = ik*It().units
    m = (m*mk.units).pm(0.1)
    return (2*np.pi*np.sqrt(1/D().value.to_base_units() *(ik.to_base_units()+It().value.to_base_units()+m.value.to_base_units()*l**2))).magnitude#, p0 = [0.1,-2.5,181]
if 1:
    T2 =  np.sum(T2, axis=0) /3  
    papt,pacov  = curve_fit(f3,l2, T2,absolute_sigma = 1, sigma = [0.006]*len(l2), p0 = [0.3,-2.5,0.01,0.04])
    popt, pcov = curve_fit(f2,l2, T2,absolute_sigma = 1, sigma = [0.006]*len(l2))#, p0 = [0.3,-2.5,0.01,0.04])
    pipt, picov = curve_fit(f4,l2, T2,absolute_sigma = 1, sigma = [0.006]*len(l2), p0 = [0.1,-2.5,181])
    
    print(SE(f2, l2, T2, popt),ChiSquare(f2, l2, T2, popt) )
    print(SE(f3, l2, T2, papt),ChiSquare(f3, l2, T2, papt) )
    print(SE(f4, l2, T2, pipt),ChiSquare(f4, l2, T2, pipt) )
    
    plt.errorbar(l2,T2, yerr = 0.006, fmt = 'g.', capsize = 3, label = 'Messwerte')   
    l2 = np.linspace(-7,1)
   # plt.plot(l2, f2(l2,*popt), label = 'Wurzelfit')
    plt.plot(l2, f3(l2,*papt),'r--', label = 'Polynomfit')
    plt.legend()

    
    
@techutil()
def g(l = l,T =T,r=r,ol=ol,ok=ok,mu=mu,m=m,phi=phi):
    return (1+2*(r/l)**2/5+ol/ok-mu/(6*m)+phi**2/8)

@techutil()
def g2(l = l,T =T):
    return 4*np.pi**2*l/T**2

@listify    
def quant(val, unit ='', err = None):
    return Q_(str(val)+unit).pm(err) if err else Q_(str(val)+unit)  

 


def error(two):
    #input list of list
    
    mean = [np.mean(a) for a in two]
    std = [np.std(a, ddof = 1) for a in two]
    stdmean = np.mean(std)
    meanerror = stdmean/np.math.sqrt(len(two[0]))
    
    return mean,std,stdmean,meanerror

def Error(ydata, xdata ):
    
    average = np.average(ydata, xdata)
    dev = [a - average for a in ydata]
    std

