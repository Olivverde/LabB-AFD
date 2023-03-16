from NFA_lab import *
import networkx as nx
import copy
# Identify available symbols
class FDA(object):

    def __init__(self) -> None:
        self.epsilon_initial_table = {}

    def table_template(self, transitions, dicc=True):
        char = ''
        elements = set()
        for i in transitions:
            elements.add(i.trans_element)
        if dicc == True:
            dicc = {elem: {} for elem in list(elements)}
        else:
            dicc = {elem: [] for elem in list(elements)}
        
        return dict(sorted(dicc.items()))

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
        self.epsilon_initial_table = copy.deepcopy(table)
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
        
                    
        # print('all_states',all_states)
        # print('afd_table',afd_table)
        
        for i in afd_table:
            for e in afd_table[i]:
                new = afd_table[i][e][0]
                afd_table[i][e] = new
        
        return all_states, afd_table

    def min_table(self, afd_table, sub_table):
        # Change values to new set values
        
        minTable = copy.deepcopy(afd_table)
        for i in afd_table:
            for k in afd_table[i]:
                for j in range(len(sub_table)):
                    if afd_table[i][k] in sub_table[j]:
                        minTable[i][k] = j
        # Finds set's mode
        index_mc = ''
        most_common = 0
        for i in sub_table[0]:
            ax = 0
            for j in sub_table[0]:
                if minTable[i] == minTable[j]:
                    ax += 1
            if ax > most_common:
                most_common = ax
                index_mc = i
        
        new_set = []
        for i in sub_table[0]:
                if minTable[i] != minTable[index_mc]:
                    new_set.append(i)
                    sub_table.insert(1,new_set)
        for i in new_set:
            sub_table[0].remove(i)
        
        return sub_table, minTable

    def minimize_afd(self):
        all_states, afd_table = self.subConstruct(transitions, acc_st)
        acceptance = len(all_states)-1
        active = True
        initial = [list(range(acceptance))]
        initial.append([acceptance])
        # print(initial)
        
        flag = 0
        while active:
            initial, minTable = self.min_table(afd_table, initial)
            if flag == len(initial):
                break
            else:
                flag = len(initial)
        
        mini = {}
        for i in range(len(initial)):
            mini[i] = {} # mini = {0:{}, 1:{}, 2:{}, 3:{}}
            for j in list(minTable[0].keys()):
                 mini[i][j] = minTable[initial[i][0]][j]  # mini = {0:{'a':x, 'b':y}, 1:{...}, 2:{...}, 3:{...}}

        return mini
                          
    def dstates_recursive(self, table, dstate_table, fp, index, sym_list):
        new_states = []
        aux_dstates = []
        dstate_table[index] = {}
        for i in sym_list:
            for j in fp:
                if table[j]['symbol'] == i: 
                    aux_dstates = list(set(aux_dstates) | set(table[j]['followpos']))
            dstate_table[index][i] = aux_dstates
            new_states.append(aux_dstates)
            aux_dstates = []
        
        return dstate_table, new_states
    
    def regular_toAFD(self, r):
        
        iteration = len(r)
        fda_stack = [char for char in r]
        index = []
        regular_table = {}
        cont = 0
        cont_nonsy = 0
        root = ''

        for i in range(iteration):
            e = fda_stack.pop(0)
            pos_template = {'nullable':None,'firstpos': [], 'lastpos':[],
                         'followpos':[], 'symbol':''}
            if e in ['.','|','*']:
                cont_nonsy += 1
                non_index = 'r'+str(cont_nonsy)
                root = non_index
                index.append(non_index)
                regular_table[non_index] = pos_template
            else:
                cont += 1
                index.append(cont)
                regular_table[cont] = pos_template
                
            if e == '.':
                pass
                e = index[i]
                c1 = index[i-2]
                c2 = index[i-1]
                c1_nullable = regular_table[c1]['nullable']
                c2_nullable = regular_table[c2]['nullable']
                
                regular_table[e]['nullable'] = c1_nullable and c2_nullable

                
                if c1_nullable:
                    regular_table[e]['firstpos'] = list(set(regular_table[c1]['firstpos'])|set(regular_table[c2]['firstpos']))
                else:
                    regular_table[e]['firstpos'] = regular_table[c1]['firstpos']
                    
                if c2_nullable:
                    regular_table[e]['lastpos'] = list(set(regular_table[c1]['lastpos'])|set(regular_table[c2]['lastpos']))
                else:
                    regular_table[e]['lastpos'] = regular_table[c2]['lastpos']

                for i in regular_table[c1]['lastpos']:
                    first_c2 = regular_table[c2]['firstpos']
                    follow_i = regular_table[i]['followpos']
                    union = list(set(follow_i)|set(first_c2))
                    regular_table[i]['followpos'] = union
                
            elif e == '|':
                
                e = index[i]
                c1 = index[i-2]
                c2 = index[i-1]
                
                regular_table[e]['nullable'] = regular_table[c1]['nullable'] or regular_table[c2]['nullable']
                regular_table[e]['firstpos'] = list(set(regular_table[c1]['firstpos'])|set(regular_table[c2]['firstpos']))
                regular_table[e]['lastpos'] = list(set(regular_table[c1]['lastpos'])|set(regular_table[c2]['lastpos']))
                                 
            elif e == '*':
                e = index[i]
                c1 = index[i-1]
                regular_table[e]['nullable'] = True
                regular_table[e]['firstpos'] = regular_table[c1]['firstpos']
                regular_table[e]['lastpos'] = regular_table[c1]['lastpos']
                
                for i in regular_table[e]['lastpos']:
                    first_n = regular_table[e]['firstpos']
                    follow_i = regular_table[i]['followpos']
                    union = list(set(follow_i)|set(first_n))
                    regular_table[i]['followpos'] = union
                
            else:  
                sym = e
                e = index[i]
                if e == 'ε':
                    regular_table[e]['nullable'] = True
                
                else:
                    regular_table[e]['nullable'] = False
                    regular_table[e]['firstpos'].append(e)
                    regular_table[e]['lastpos'].append(e)
                    regular_table[e]['symbol'] = sym  
        
        sym_set = set()
        for inner_dicc in regular_table.values():
            symbol = inner_dicc.get("symbol")
            if symbol and symbol.strip():
                sym_set.add(symbol.strip())
        sym_list = sorted(list(sym_set))
        sym_list.remove('#')
        
        dstates = [regular_table[root]['firstpos']]
        dstates_table = {}
        new_states = []
        for fp in dstates:
            fp_ind = dstates.index(fp)
            dstates_table, new_states = self.dstates_recursive(regular_table, dstates_table, fp, fp_ind, sym_list)
            
            for ns in new_states:
                if not ns in dstates:
                    dstates.append(ns)
        
        for i in sym_list:
            for fp in dstates_table:
                stt = dstates_table[fp][i]
                dstates_table[fp][i] = dstates.index(stt)
                
        
        return dstates_table

    def graph(self, dicc):
        
        G = nx.DiGraph()
        
        for i in dicc:
            for j in dicc[i]:
                G.add_node(i)
                G.add_node(dicc[i][j])
                # Agregar una transición del nodo 0 al nodo 1 con el símbolo 'a'
                G.add_edge(i, dicc[i][j], label=j)
        
        pos = nx.spring_layout(G)
        nx.draw_networkx_nodes(G, pos)
        nx.draw_networkx_edges(G, pos)
        nx.draw_networkx_edge_labels(G, pos)
        nx.draw_networkx_labels(G, pos)
        plt.show()
            
    def afd_simulation(self, w, dicc):
        
        initial_state = 0
        s = initial_state
        final_state = len(dicc)-1
        for c in w:
            s = dicc[s][c]
        if s == final_state:
            return 'PASS'
        else:
            return 'FAIL'

    def afn_simulation(self, w):
        
        epsilon_table = self.epsilon_initial_table
        F = len(epsilon_table)-1
        s = epsilon_table[0]['ε']
        ax = []
        ax2 = []
        
        for c in w:
            for i in s:
                e = epsilon_table[i][c]
                if isinstance(e,int):
                    ax.append(epsilon_table[i][c])
            for j in ax:
                ax2 = list(set(ax2) | set(epsilon_table[j]['ε']))
            s = ax2
            ax = []
            ax2 = []
        if F in s:
            return 'PASS'
        else:
            return 'FAIL' 

r = ['(a|b)*abb']
q = 0
re_list = {'basic':r[q], 'regular':r[q]+'#'}
re = re_list['basic']

lib = Libs(re)
postfix = lib.get_postfix()
print('---------------------------')
print('TRADUCCION:',lib.get_printable_trans())
print('POSTFIX:',lib.get_printable_postfix())
print('---------------------------')

nfa = NFA()
fda = FDA()

postfix_fda = 'ab|*a.b.b.#.'

w = 'aabb'
nfa.thompson(postfix)

# states, nodes, transitions, etc. NFA in a nutshell
transitions = nfa.get_transitions()
acc_st = nfa.get_acceptance_state()+1
a,b = fda.subConstruct(transitions, acc_st)
print('NFA')
print('-----------------------------------------------------------')
print('NFA set →',nfa.get_standarized_trans())
print('W =',w,' → Simulation Status:',fda.afn_simulation(w))
print('-----------------------------------------------------------')
print('FDA')
print('-----------------------------------------------------------')
print('Subconstruction DFA Set →',b)
print('W =',w,' → Simulation Status:',fda.afd_simulation(w,b))
print('-----------------------------------------------------------')
fda.graph(b)
print('Subconstruction DFA Set Minimized →',fda.minimize_afd())
print('W =',w,' → Simulation Status:',fda.afd_simulation(w,fda.minimize_afd()))
print('-----------------------------------------------------------')
fda.graph(fda.minimize_afd())
print('Direct FDA Construction →',fda.regular_toAFD(postfix_fda))
print('W =',w,' → Simulation Status:',fda.afd_simulation(w,fda.regular_toAFD(postfix_fda)))
print('-----------------------------------------------------------')
fda.graph(fda.regular_toAFD(postfix_fda))
