from de_bruijn import DeBruijnGraph
from k_mers import frequency_error_correction, KMerIterator
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

    # Todas las lecturas
    raw_reads = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        for k_mer in KMerIterator(read, k):
            raw_reads.append(k_mer)
    write_assembly(raw_reads, 'raw.fasta', 'raw reads')

    # Lecturas con q-score >= 30 y longitud entre 100 y 150
    quality_reads = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        if min(quality) >= '?' and 100 <= len(read) <= 150:
            for k_mer in KMerIterator(read, k):
                quality_reads.append(k_mer)
    write_assembly(quality_reads, 'quality.fasta', 'quality reads')

    # Lecturas con q-score >-0 30 y longitud entre 100 y 150
    # Solo se toman los k-meros que aparecen más de 4 veces en las lecturas
    frequency_reads = []
    for read, quality in FastqIterator('bacata_reads_1.fastq'):
        if min(quality) >= '?' and 100 <= len(read) <= 150:
            for k_mer in KMerIterator(read, k):
                frequency_reads.append(k_mer)
    frequency_reads = frequency_error_correction(frequency_reads, 4)
    write_assembly(frequency_reads, 'frequency.fasta', 'frequency reads')

    # Lecturas simuladas
    genome = read_fasta('bacata.fasta')
    simulated_reads = []
    for read in simulate_sequencing(genome, 54495, 100, 150, linear=False):
        for k_mer in KMerIterator(read, k):
            simulated_reads.append(k_mer)
    write_assembly(simulated_reads, 'simulated.fasta', 'simulated reads')
