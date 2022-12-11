def counts_dict(k_mers):
    '''
    Regresa un diccionario que asocia a cada k-mero\n
    la cantidad de veces que aparece en k_mers
    '''
    counts = {}
    for k_mer in k_mers:
        counts[k_mer] = counts.get(k_mer, 0) + 1
    return counts


def frequency_error_correction(k_mers, threshold):
    '''
    Recibe una lista de k-meros y regresa una nueva lista\n
    que contiene a los k-meros que aparecen más de threshold veces
    '''
    corrected = []
    counts = counts_dict(k_mers)
    while len(k_mers) > 0:
        k_mer = k_mers.pop()
        if counts[k_mer] >= threshold:
            corrected.append(k_mer)
    corrected.reverse()
    return corrected


def transition_transversion_mutations(n):
    '''
    Recibe una base t regresa una lista con las\n
    posibles transiciones y transversiones de esa base
    '''
    if n == 'A':
        return ['C', 'G', 'T']
    elif n == 'C':
        return ['A', 'G', 'T']
    elif n == 'G':
        return ['A', 'C', 'T']
    return ['A', 'C', 'G']


def hamming_distance_corrected(k_mer, counts, threshold):
    '''
    Recibe un k-mero y regresa el k-mero corregido\n
    usando la distancia de Hamming
    '''
    max = k_mer
    max_count = counts[k_mer]
    for index in range(len(k_mer)):
        for base in transition_transversion_mutations(k_mer[index]):
            neighbor = k_mer[:index] + base + k_mer[index + 1:]
            neighbor_count = counts.get(neighbor, 0)
            if neighbor_count >= threshold:
                return neighbor
            elif neighbor_count > max_count:
                max = neighbor
    return max


def hamming_distance_error_correction(k_mers, threshold):
    '''
    Recibe una lista de k-meros y regresa una nueva lista\n
    que contiene a los k-meros corregidos usando la distancia de Hamming
    '''
    corrected = []
    counts = counts_dict(k_mers)
    while len(k_mers) > 0:
        k_mer = k_mers.pop()
        if counts[k_mer] < threshold:
            k_mer = hamming_distance_corrected(k_mer, counts, threshold)
        corrected.append(k_mer)
    corrected.reverse()
    return corrected
