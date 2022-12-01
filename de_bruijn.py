class DeBruijnNode:
    '''
    Clase que representa un nodo de una gráfica de De Bruijn
    '''

    def __init__(self):
        '''
        Inicializa un nuevo nodo sin aristas que entran ni salen
        '''
        self.edges = {}
        self.balance = 0


class DeBruijnGraph:
    '''
    Clase que representa una gráfica de De Bruijn para ensamblar un genoma\n
    a partir de los k-meros de lecturas de secuenciación
    '''

    def __init__(self, k_mers):
        '''
        Inicializa la gráfica a partir de los k-meros especificados
        '''
        self.nodes = {}
        for k_mer in k_mers:
            prefix = k_mer[:-1]
            suffix = k_mer[1:]
            self.add_node(prefix)
            self.add_node(suffix)
            self.add_edge(prefix, suffix)

    def add_node(self, string):
        '''
        Agrega un nodo con la cadena especificada
        '''
        self.nodes.setdefault(string, DeBruijnNode())

    def add_edge(self, prefix, suffix):
        '''
        Agrega una arista dirigida del nodo con la cadena prefix al nodo con la cadena suffix
        '''
        node = self.nodes[prefix]
        node.balance += 1
        node.edges[suffix] = node.edges.get(suffix, 0) + 1
        node = self.nodes[suffix]
        node.balance -= 1

    def has_eulerian_path(self):
        '''
        Regresa una tupla donde la primera entrada es True si la gráfica tiene un camino euleriano y la segunda entrada es la cadena del nodo donde empieza tal ciclo
        '''
        start = None
        end = None
        for string, node in self.nodes.items():
            if node.balance == 1 and not start:
                start = string
            elif node.balance == -1 and not end:
                end = string
            elif node.balance != 0:
                return False, None
        if not (start and end):
            return False, None
        return True, start

    def has_eulerian_cycle(self):
        '''
        Regresa True si la gráfica tiene un ciclo euleriano
        '''
        for node in self.nodes.values():
            if node.balance != 0:
                return False
        return True

    def continue_to_next_node(self, node, path, stack):
        '''
        Se manda a llamar durante el algoritmo de Hierholzer en el caso\n
        en el que todavía hay aristas por explarar que salen del nodo actual
        '''
        string, weight = node.edges.popitem()
        path.append(string)
        if weight > 1:
            node.edges[string] = weight - 1
        elif len(node.edges) > 0:
            stack.append(node)

    def backtrack(self, path, stack, genome):
        '''
        Se manda a llamar durante el algoritmo de Hierholzer en el caso\n
        en el que ya no hay más aristas que salgan del nodo actual
        '''
        node = stack.pop() if len(stack) > 0 else None
        while len(path) > 0:
            string = path[-1]
            if self.nodes[string] is node:
                return
            else:
                path.pop()
                genome.append(string[0])

    def assembly(self, string):
        '''
        Ensambla el genoma construyendo un camino euleriano con el algoritmo de Hierholzer\n
        El camino empieza en el nodo con la cadena especificada\n
        Este método elimina las aristas de la gráfica
        '''
        path = [string]
        stack = []
        genome = []
        while len(path) > 0:
            string = path[-1]
            node = self.nodes[string]
            if len(node.edges) > 0:
                self.continue_to_next_node(node, path, stack)
            else:
                self.backtrack(path, stack, genome)
        genome.reverse()
        return ''.join(genome)
