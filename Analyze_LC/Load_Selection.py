import numpy as np
import operator


def get_operator_fn(op):
    return {
        '+' : operator.add,
        '-' : operator.sub,
        '*' : operator.mul,
        '/' : operator.div,
        '%' : operator.mod,
        '^' : operator.xor,
        '>=': operator.ge,
        '>': operator.gt,
        '<=': operator.le,
        '<' : operator.lt,
        '==': operator.eq,
        }[op]



test=np.loadtxt('Sela_with_colorcut.txt',dtype={'names': ('cut','comp','val','type'),'formats': ('S15','S2','S8','S5')})

print(test)

sel=[]
for cut in test:
    thecut=[]
    for val in cut['cut'].split('+'):
        thecut.append(val)
    print(cut['type'])
    if cut['type'] != 'str':
        sel.append((thecut,get_operator_fn(cut['comp']),eval(cut['type']+'('+cut['val']+')')))
    else:
        sel.append((thecut,get_operator_fn(cut['comp']),cut['val']))

print(sel)
