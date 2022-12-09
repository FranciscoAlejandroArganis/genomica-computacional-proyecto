def get_k_mer(sequence, k, start):
    '''
    Regresa el k-mero de la cadena sequence que empieza en el índice start
    '''
    end = start + k
    if end <= len(sequence):
        return sequence[start:end]
    return sequence[start:] + sequence[:end - len(sequence)]


def frequency_error_correction(k_mers, threshold):
    '''
    Recibe una lista de k-meros y regresa una nueva lista\n
    que contiene a los k-meros que aparecen más de threshold veces
    '''
    res = []
    freq = {}
    total = 0
    for k_mer in k_mers:
        freq[k_mer] = freq.get(k_mer, 0) + 1
        total += 1
    while len(k_mers) > 0:
        k_mer = k_mers.pop()
        if freq[k_mer] >= threshold:
            res.append(k_mer)
    return res


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
