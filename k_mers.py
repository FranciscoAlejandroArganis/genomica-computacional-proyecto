class KMerIterator:
    '''
    Clase que permite iterar por todos los k-meros de una cadena
    '''

    def __init__(self, string, k):
        '''
        Inicializa un nuevo iterador de k-meros para la cadena especificada
        '''
        self.string = string
        self.k = k

    def __iter__(self):
        '''
        Regresa el iterador
        '''
        self.index = 0
        return self

    def __next__(self):
        '''
        Regresa el siguiente k-mero durante la iteración
        '''
        if self.index <= len(self.string) - self.k:
            k_mer = self.string[self.index:self.index + self.k]
            self.index += 1
            return k_mer
        raise StopIteration
