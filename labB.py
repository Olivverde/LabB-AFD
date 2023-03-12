from NFA_lab import *
re_list = ['(a|b)*abb']
re = re_list[0]

lib = Libs(re)
postfix = lib.get_postfix()
# print('---------------------------')
# print('TRADUCCION:',lib.get_printable_trans())
# print('POSTFIX:',lib.get_printable_postfix())
# print('---------------------------')

nfa = NFA()
nfa.thompson(postfix)
# states, nodes, transitions, etc. NFA in a nutshell
transitions = nfa.get_transitions()
acc_st = nfa.get_acceptance_state()+1

# Identify available symbols
def table_template(transitions, dicc=True):
    char = ''
    elements = set()
    for i in transitions:
        elements.add(i.trans_element)
    if dicc == True:
        dicc = {elem: {} for elem in list(elements)}
    else:
        dicc = {elem: [] for elem in list(elements)}
    
    return dicc

def epsilon_recursive(init_origin, i):
    for j in transitions.sub_transitions:
        
        init_dg = j.local_init_node.id
        end_dg = j.local_end_node.id
        elem_dg = j.trans_element
        if init_dg == i:
            if elem_dg != 'ε' and init_origin == None:
                table[i][elem_dg] = end_dg
            elif elem_dg == 'ε':
                try:
                    table[init_origin][elem_dg].append(end_dg)
                except:
                    table[i][elem_dg].append(end_dg)
                    init_origin = i
                epsilon_recursive(init_origin, end_dg)
                

def state_recursive(table, states_list, initial_list, element=None, aux_list=[]):
    # Check if there is a current recursivity cycle 
    if element == None:
        ele_aux = states_list
    elif element != None:
        ele_aux = [element]
    for k in ele_aux:
        for j in initial_list:
            if isinstance(table[j][k], int) and element == None:
                    aux_list = state_recursive(table, states_list, [table[j][k]], k, aux_list)
            elif isinstance(table[j][k], int) and element != None:
                union = set(aux_list) | set(table[j]['ε'])
                aux_list = list(union)
                return aux_list         
            elif (not isinstance(table[j][k], int)) and (element != None):
                union = set(aux_list) | set(table[j]['ε'])
                aux_list = list(union)
                return aux_list
            
        if not aux_list in states_list[k]:
            states_list[k].append(aux_list)
        aux_list = []
    return states_list

# Transition Table Template
table = table_template(transitions.sub_transitions)
table = {i: table.copy() for i in range(acc_st)}
for state in table.copy():
    table[state]['ε'] = [state]

# Generate AFD transitions table
for i in range(acc_st):
    
    epsilon_recursive(None, i)

# Generate AFD states table
states_list = table_template(transitions.sub_transitions, False)
states_list.pop('ε')
initial_list = table[0]['ε']

all_states = [initial_list]
for i in all_states:
    states_list = (state_recursive(table, states_list, i))
    for k in states_list:
        for j in states_list[k]:
            if not j in all_states:
                all_states.append(j)
print('all_states',all_states)

# print(table)
            
        