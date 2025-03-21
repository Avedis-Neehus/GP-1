
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 19:20:35 2018
@author: Avedis Neehus
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#from statsmodels.stats.weightstats import DescrStatsW as data
#from uncertainties import ufloat as u, umath
from pint import  UnitRegistry  
import sympy as sym
import inspect

ur = UnitRegistry()
Q_ = ur.Quantity
Q_.pm = Q_.plus_minus
Q_.mag= Q_.magnitude 

class attriList(list):
    
    ''' Custom list defined to have acces to attributes of elements of lists without having to expicitly lopp through the list
    everytime. Doesnt work for methods for somw reason'''
    
    def __getattr__(self,attr):
        
            return attriList(item.__getattribute__(attr) for item in self)
        
def listify(f):
    '''Decorator to shorten map syntax. Only vectorizes *args'''    
    def wrapping(*args,**kwargs):        
        return attriList(map(lambda *arg : f(*arg,**kwargs), *args))
    
    return wrapping
@listify    
def quant(val, unit ='', err = None):
    return Q_(str(val)+unit).pm(err) if err else Q_(str(val)+unit)
   
class mpi:
    '''Switches behavior between symp and numpy functions/constants, which is needed for correct behavior of errify and
    texify. Instantiate with sympy name of function/constant as string and if the name is different in numpy use it as 
    a second argument'''
    
    symode = True
    
    def __init__(self, symb, nump = None):
        if type(symb) ==str:
            self.sym = eval('sym.'+symb)
            self.np  =  eval('np.'+nump) if nump else eval('np.' + symb)
        elif isinstance(symb,(float,int)):
            self.sym = sym.Symbol(str(symb))
            self.np = symb            
        else:
            self.sym = symb
            self.np = float(symb)
            
    def __repr__(self):
        return str(self.sym)    
                
    def __call__(self,calle):
        
        if self.symode:
            return (self.sym(calle.sym)) if hasattr(calle,'sym') else (self.sym(calle))
        else: return (self.np(calle.np)) if hasattr(calle,'sym') else (self.np(calle))
    
    def  __mul__(self, other):
       
       if self.symode:
           return (self.sym * other)
       else: 
           return (self.np*other)
       
    def  __truediv__(self, other):
       
       if self.symode:
           return self.sym / other
       else: 
           return self.np/other
       
    def __add__(self, other):
       
       if self.symode:
           return self.sym + other
       else: 
           return self.np+other
    def __sub__(self, other):
        
        if self.symode:
           return self.sym - other
        else: 
           return self.np-other 
    def __pow__(self,other):
        
        if self.symode:
           return self.sym** other
        else: 
           return self.np**other
       
pi = mpi('pi')
sin = mpi('sin')  
cos = mpi('cos') 
tan = mpi('tan')  
        
def texify( *varnames,fname = None, **cnames): 
    '''Decorator to print latex code for a function. The function must be pure and consist 
    of elementary operations. The variables should all be pint objects. Constants have to be provided as keywordarguments, they are static.'''
    def deco(f):        
        def function(*args):
            mpi.symode = 1
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
            
            mpi.symode = 0
            print('\\begin{align*}')
            print(entity + ' &= '+ sym.latex(symexpr) + r'\\')
            
            print('&=' + sym.latex(numexpr, mul_symbol = 'dot') + '=', end = '')
           # try:
            result = f(*args)
            print( '{:Lx.02f}'.format(result))
            #except: pass
            print('\\end{align*}')
            
        return function
    return deco

def errify( *varnames, fname = None, **cnames):
    
    def deco(f):
        
        def function(*args):
            
            mpi.symode = 1
            
            entity  = fname if fname else f.__name__
            
            #get the symbols for arguments in func definition
            argspec =  inspect.getfullargspec(f)            
            def_args = argspec.args[:-len(argspec.defaults)] if argspec.defaults else argspec.args
            def_kwargs = argspec.args[-len(argspec.defaults):] if argspec.defaults else {}
            
            #Create sympy symbols
            symargs   = tuple(sym.Symbol(varname)  for varname in varnames)  if varnames else tuple(sym.symbols(def_args))
            symkwargs = tuple(sym.Symbol(cname) for cname in cnames.values())if cnames   else tuple(sym.symbols(def_kwargs))
            symall = symargs+symkwargs
            
            allarg = args+f.__defaults__ if f.__defaults__ else args
            
            #Get symbols with attached uncertainty (Measurements)
            symprecise = [symbol for symbol,arg in zip(symall,allarg) if hasattr(arg, 's')]
            
            #Calls function with sympy tokens and returns sympy expression
            expr = f(*symall)#, **symkwargs)            
            derivatives = [expr.diff(symbol) for symbol in symprecise]
            
            #Same as above but with numerical symbols
            inprecise = [sym.Symbol('{:Lx}'.format(arg.error)) for arg in allarg if hasattr(arg, 's')]
            numargs = tuple(sym.Symbol('{:Lx}'.format(arg.value)) if  hasattr(arg, 's') else sym.Symbol('{:Lx}'.format(arg)) for arg in allarg  )
            numexpr = f(*numargs)
            numderivs =  [numexpr.diff(symbol) for symbol,val in zip(numargs,allarg) if hasattr(val, 's')]
            
            #print the result
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
            
             
            
            mpi.symode = 0       
            value = f(*args)
            print(' &= ' , '{:Lx.03f}'.format(value.s * value.units))  
        
            print(r'\end{align*}')
            
        return function 
    return   deco         
        

def techutil(*varnames,fname = None, **cnames):
    '''Decorator which attaches the texified and errified versions of the function as attributes to the function'''
    def deco(f):
                
        def function(*args):            
            return f(*args)    
        
        function.gauß = errify(*varnames,fname = fname, **cnames)(f)
        function.tex = texify(*varnames,fname= fname , **cnames)(f)   
        
        return function    
    return deco


def tabelize(data, collab = None, rowlab = None, hline = 0):
    '''Takes a matrix and spits out a the unstylized code for a latex table.
    collab and rowlab are iterables with the names of the column, row'''
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
  



def barplot(data,xaxis, xlabel = None, ylabel = None, maxy= None, labels = None):
     
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

def mprint(q):
    print( '{:Lx}'.format(q))

#some function for regression error estimation
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

if 0:
    d = Q_('1mm')
    
    Bp1 = Q_('1200mT').pm(5) #15V Abstandfehler 50
    gepI1=quant([35,22,8,49,15,29,43], unit = 'mA') #mA
    gepU1= quant([26,21,16,30,18,23,28], unit = 'V')#*0.1 #V
    
    Bp2 = Q_('753mT').pm(5)#10V
    gepI2= quant([8,15,22,29,36,43,49], unit = 'mA') #mA
    gepU2 =quant([13,15,18,19,21,23,24], unit = 'V')#*0.1 #V 
    
    Bn1 = Q_('1185mT').pm(5)#15V
    GenI1 = quant([1,8,15,22,29,36,43,50], unit = 'mA')#mA
    GenU1 = quant([43,51,60,68,76,85,92,100], unit = 'V')#*0.1#V
    
    Bn2 = Q_('815mT').pm(5)#10
    GenI2 = quant([1.5,8,16,22,31,36,43,49], unit = 'mA')#mA
    GenU2 = quant([34,40,47,53,61,65,71,78], unit = 'V') #V
    
    @techutil()
    def Rh(B,I,U, d = d):
        return (U*d)/(B*I)
    
    print(Rh.gauß(Bp1,gepI1[0],gepU1[0]*0.1))

    
            
