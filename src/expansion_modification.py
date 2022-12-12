from error_correction import transition_transversion_mutations
from random import random, randint


def exp_mod_itr(sequence, p):
    '''
    Realiza la siguiente iteración en el proceso de expansión-modificación\n
    Con probabilidad p ocurre una mutación en la secuencia\n
    Con probabilidad 1 - p se repiten las últimas 100 bases en la secuencia
    '''
    if random() <= p:
        index = randint(0, len(sequence) - 1)
        sequence[index] = transition_transversion_mutations(sequence[index])[
            randint(0, 2)]
    else:
        length = len(sequence)
        index = length - 100
        if index < 0:
            index = 0
        while index < length:
            sequence.append(sequence[index])
            index += 1


def synthetic_genome(iterations, p):
    '''
    Regresa un genoma sintético generado por un\n
    proceso de expansión-modificación con probabilidad de mutación p
    '''
    sequence = [['A', 'C', 'G', 'T'][randint(0, 3)]]
    for _ in range(iterations):
        exp_mod_itr(sequence, p)
    return ''.join(sequence)
