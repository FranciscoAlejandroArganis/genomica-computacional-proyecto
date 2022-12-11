def get_k_mer(sequence, k, start):
    '''
    Regresa el k-mero de la cadena sequence que empieza en el índice start
    '''
    end = start + k
    if end <= len(sequence):
        return sequence[start:end]
    return sequence[start:] + sequence[:end - len(sequence)]


class KMerIterator:
    '''
    Clase que permite iterar por todos los k-meros de una cadena\n
    Si linear es False, se considera que la cadena es circular
    '''

    def __init__(self, string, k, linear=True):
        '''
        Inicializa un nuevo iterador de k-meros para la cadena especificada
        '''
        self.string = string
        self.k = k
        self.limit = len(string) - k if linear else len(self.string) - 1

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
        if self.index <= self.limit:
            k_mer = get_k_mer(self.string, self.k, self.index)
            self.index += 1
            return k_mer
        raise StopIteration
