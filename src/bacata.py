from de_bruijn import DeBruijnGraph
from k_mers import KMerIterator
from error_correction import frequency_error_correction, hamming_distance_error_correction
from fasta import read_fasta, write_fasta, simulate_sequencing, FastqIterator


def write_assembly(k_mers, path, header, unique=True):
    '''
    Escribe en un archivo fasta el resultado de ensamblar el genoma con los k-meros especificados\n
    Si unique es True también escribe en un segundo archivo el ensamble con el conjunto de k-meros únicos
    '''
    de_bruijn_graph = DeBruijnGraph(k_mers)
    for string in de_bruijn_graph.nodes:
        break
    assembled_genome = de_bruijn_graph.assembly(string)
    write_fasta(path, header, assembled_genome)
    if unique:
        de_bruijn_graph = DeBruijnGraph(set(k_mers))
        for string in de_bruijn_graph.nodes:
            break
        assembled_genome = de_bruijn_graph.assembly(string)
        write_fasta('u-' + path, 'unique ' + header, assembled_genome)


# Script que genera ensambles del genoma del fago Bacata de Xylella
# a partir de los datos de secuenciación reales y simulados
if __name__ == '__main__':
    k = 35
    genome = read_fasta('bacata.fasta')

    # Todas las lecturas
    k_mers = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        for k_mer in KMerIterator(read, k):
            k_mers.append(k_mer)
    write_assembly(k_mers, 'raw.fasta', 'raw reads')

    # Lecturas con buena calidad
    k_mers = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        if min(quality) >= '?' and 100 <= len(read) <= 150:
            for k_mer in KMerIterator(read, k):
                k_mers.append(k_mer)
    write_assembly(k_mers, 'quality.fasta', 'quality reads')

    # Lecturas corregidas con frecuencias
    k_mers = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        if min(quality) >= '?' and 100 <= len(read) <= 150:
            for k_mer in KMerIterator(read, k):
                k_mers.append(k_mer)
    k_mers = frequency_error_correction(k_mers, 4)
    write_assembly(k_mers, 'frequency.fasta', 'frequency reads')

    # Lecturas corregidas con distancia de Hamming
    k_mers = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        if min(quality) >= '?' and 100 <= len(read) <= 150:
            for k_mer in KMerIterator(read, k):
                k_mers.append(k_mer)
    k_mers = hamming_distance_error_correction(k_mers, 4)
    write_assembly(k_mers, 'hamming.fasta', 'Hamming distance reads')

    # Lecturas simuladas
    k_mers = []
    for read in simulate_sequencing(genome, 54495, 100, 150, linear=False):
        for k_mer in KMerIterator(read, k):
            k_mers.append(k_mer)
    write_assembly(k_mers, 'simulated.fasta', 'simulated reads')
