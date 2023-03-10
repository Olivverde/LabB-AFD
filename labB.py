import os
from LAB_A.labA import *

re_list = ['(a|b)*abb']
re = re_list[0]

lib = Libs(re)
postfix = lib.get_postfix()
print('---------------------------')
print('TRADUCCION:',lib.get_printable_trans())
print('POSTFIX:',lib.get_printable_postfix())
print('---------------------------')
nfa = NFA()
nfa.thompson(postfix)

