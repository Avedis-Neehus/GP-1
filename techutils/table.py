# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 19:20:35 2018

@author: Elios
"""

import numpy as np
import matplotlib.pyplot as plt
#from statsmodels.stats.weightstats import DescrStatsW as data
from uncertainties import ufloat as u
from pint import  UnitRegistry  
import sympy as sym

ur = UnitRegistry()
Q_ = ur.Quantity
Q_.pm = Q_.plus_minus

def listify(f):
        
    def wrapping(*args,**kwargs):        
        return list(map(lambda *arg : f(*arg,**kwargs), *args))
    
    return wrapping


def errify(entity, *varnames, **cnames):
    
    def deco(f):
        
        def function(*args):
            
            symargs   = tuple([sym.Symbol(varname) for varname in varnames])
            symkwargs = dict((key,sym.Symbol(val)) for key,val in cnames.items())
            symall = symargs+tuple(symkwargs.values())
            allarg = args+f.__defaults__
            
            print(symall,allarg)
           # argindex = [ i for i in range(args) if hasattr(args[i], 'std_dev')]
            symprecise = [symbol for symbol,arg in zip(symall,allarg) if hasattr(arg, 's')]
            inprecise = [arg for arg in allarg if hasattr(arg, 's')]
            #syment    = sym.symbols(entity)#(*inprecise) 
           
            
            expr = f(*symargs, **symkwargs)
            value = f(*args)
            derivatives = [expr.diff(symbol) for symbol in symprecise]
            
           # evalderiv   = [sym.lambdify(*symall, deriv) for deriv in derivatives]
            newargs = tuple('{:Lx.02f}'.format(arg.value) if hasattr(arg, 's') else '{:Lx.02f}'.format(arg) for arg in allarg )
            
            
            
            print('Der Fehler von ' + entity+ ' ist:')
            print('\\begin{align*}')
            
            print('\\Delta ' + entity + str(symall) + '&= \\sqrt{', end =  '')
            
            for symbol in symprecise:
                print('\\left(' + '\\frac{\\partial ' + entity + '}{\\partial ' + sym.latex(symbol) + '}'  + '\\Delta ' + sym.latex(symbol) +  ' \\right)^{2}' + '+', end = ' ' ) 
            
                
            print(r' } \\ &= ')
            
            print(r'\sqrt{', end = '')
            for val,derivativ in zip(inprecise,derivatives):
                print('\\left(' + sym.latex(derivativ) + ' \\cdot ' + '{:Lx.02f}'.format(val.error)+'\\right)^{2} +', end = ' ' )
            print(r'} \\')    
            
            print(entity + '\\left(')
            
            for arg in newargs :
                print(arg, end = ',')
            print(' \\right)')# '& =' )   
                      
            #print(r'\sqrt{', end = '')
            
            #for  val,error in zip(inprecise,value.error_components().values()) :
                                  
            #    print('\\left(' + '{:Lx.03f}'.format((error*value.units)/(val.error)) + ' \\cdot' + '{:Lx.02f}'.format(val.error) + '\\right)^{2} + ')
            #print(r'}\\')  
                
            print(' &= ' , '{:Lx.03f}'.format(value.s * value.units))   
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
    
@errify('f', s = 's', e ='e')
def brenn(    
    s = Q_('44.2cm').pm(0.2),
    e = Q_('26.5cm').pm(0.2) ):
    
    f = (s**2-e**2)/(4*s)
    return f
    
brenn()


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

