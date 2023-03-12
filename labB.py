from NFA_lab import *

# Identify available symbols
class FDA(object):

    def table_template(self, transitions, dicc=True):
        char = ''
        elements = set()
        for i in transitions:
            elements.add(i.trans_element)
        if dicc == True:
            dicc = {elem: {} for elem in list(elements)}
        else:
            dicc = {elem: [] for elem in list(elements)}
        
        return dicc

    def epsilon_recursive(self, table, init_origin, i):
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
                    self.epsilon_recursive(table, init_origin, end_dg)
                    
    def state_recursive(self, table, states_list, initial_list, element=None, aux_list=[]):
        # Check if there is a current recursivity cycle 
        if element == None:
            ele_aux = states_list
        elif element != None:
            ele_aux = [element]
        for k in ele_aux:
            for j in initial_list:
                if isinstance(table[j][k], int) and element == None:
                        aux_list = self.state_recursive(table, states_list, [table[j][k]], k, aux_list)
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

    def subConstruct(self, transitions, acc_st):
    # Transition Table Template
        table = self.table_template(transitions.sub_transitions)
        table = {i: table.copy() for i in range(acc_st)}
        for state in table.copy():
            table[state]['ε'] = [state]

        # Generate AFD transitions table
        for i in range(acc_st):
            self.epsilon_recursive(table, None, i)

        # Generate AFD states table
        states_list = self.table_template(transitions.sub_transitions, False)
        states_list.pop('ε')
        initial_list = table[0]['ε']

        afd_table = {}
        trans_table = {}
        all_states = [initial_list]
        states_cont = 0
        for i in all_states:
            trans_table[states_cont] = i
            states_list = (self.state_recursive(table, states_list, i))
            for k in states_list:
                for j in states_list[k]:
                    if not j in all_states:
                        all_states.append(j)
            afd_table[states_cont] = states_list
            states_list = self.table_template(transitions.sub_transitions, False)
            states_list.pop('ε')
            states_cont += 1


        for i in afd_table:
            for e in afd_table[i]:
                for t in trans_table:
                    if trans_table[t] == afd_table[i][e][0]:
                        afd_table[i][e][0] = t
                    
        print('all_states',all_states)
        print('afd_table',afd_table)

re_list = ['(a|b)*abb']
re = re_list[0]

lib = Libs(re)
postfix = lib.get_postfix()
# print('---------------------------')
# print('TRADUCCION:',lib.get_printable_trans())
# print('POSTFIX:',lib.get_printable_postfix())
# print('---------------------------')

nfa = NFA()
fda = FDA()
nfa.thompson(postfix)
# states, nodes, transitions, etc. NFA in a nutshell
transitions = nfa.get_transitions()
acc_st = nfa.get_acceptance_state()+1

fda.subConstruct(transitions, acc_st)
            
        