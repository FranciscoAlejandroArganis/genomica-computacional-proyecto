from random import randrange
from k_mers import get_k_mer


def read_fasta(path):
    '''
    Regresa una cadena con la secuencia leida de un archivo fasta\n
    Se asume que el archivo solo contiene una secuencia
    '''
    with open(path, 'r') as file:
        lines = []
        file.readline()
        while True:
            line = file.readline()
            if not line:
                break
            lines.append(line[:-1])
    return ''.join(lines)


def write_fasta(path, header, sequence):
    '''
    Escribe en un archivo fasta la secuencia especificada, con la cadena header en la primera línea
    '''
    with open(path, 'a') as file:
        file.write('> ')
        file.write(header)
        file.write('\n')
        for fragment in [sequence[index:index + 70] for index in range(0, len(sequence), 70)]:
            file.write(fragment + '\n')


def simulate_sequencing(genome, n, min_len, max_len, linear=True):
    '''
    Simula la secuenciación de un genoma\n
    Regresa n lecturas con tamaño entre min_len y max_len\n
    '''
    limit = len(genome) - length if linear else len(genome)
    reads = []
    for _ in range(n):
        length = randrange(min_len, max_len)
        start = randrange(0, limit)
        reads.append(get_k_mer(genome, length, start))
    return reads


class FastqIterator:
    '''
    Clase que permite iterar por todas las lecturas en un archivo fastq
    '''

    def __init__(self, path):
        self.path = path

    def __iter__(self):
        '''
        Regresa el iterador
        '''
        self.file = open(self.path, 'r')
        return self

    def __next__(self):
        '''
        Regresa una tupla donde la primer entrada es la siguiente lectura del archivo y la segunda es la calidad
        '''
        line = self.file.readline()
        if line:
            line = self.file.readline()
            self.file.readline()
            quality = self.file.readline()
            return line[:-1], quality[:-1]
        self.file.close()
        raise StopIteration
