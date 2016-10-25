from __init__ import *

objt = tdpy.util.gdatstrt()
objt.varb = 1.

arry = zeros(1000)

def funcarry():

    temp = arry[0]


def funcobjt():

    temp = objt.varb


tdpy.util.time_func_verb(funcarry)
tdpy.util.time_func_verb(funcobjt)


