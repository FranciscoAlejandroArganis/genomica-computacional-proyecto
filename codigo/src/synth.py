from de_bruijn import DeBruijnGraph
from k_mers import KMerIterator
from expansion_modification import synthetic_genome
from fasta import read_fasta, write_fasta, simulate_sequencing

# Script que genera ensambles del genoma sintético
# a partir de datos de secuenciación simulados
# con diferentes valores de k
if __name__ == '__main__':
    genome = synthetic_genome(1000, .75)
    write_fasta('synth.fasta', 'synthetic genome', genome)
    genome = read_fasta('synth.fasta')

    # Ensamble del genoma sintético con k desde 20 hasta 50
    reads = simulate_sequencing(genome, 50000, 100, 150)
    for k in range(20, 51, 5):
        k_mers = []
        for read in reads:
            for k_mer in KMerIterator(read, k):
                k_mers.append(k_mer)
        de_bruijn_graph = DeBruijnGraph(k_mers)
        assembled_genome = de_bruijn_graph.assembly(genome[:k - 1])
        write_fasta('asm_synth_' + str(k) + '.fasta',
                    'assembled synthetic genome k = ' + str(k), assembled_genome)
