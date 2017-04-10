import argparse
import itertools
from tqdm import tqdm


class Edge:
    show_sequences = False

    def __init__(self, v1, v2):
        assert isinstance(v1, Vertex) and isinstance(v2, Vertex)
        self.v1 = v1
        self.v2 = v2
        self.seq = v1.seq + v2.seq[-1]
        self.cov = 1

    def inc_coverage(self):
        self.cov += 1

    def __len__(self):
        return len(self.seq)

    def merge(self, following_edge):
        self.v2 = following_edge.v2
        self.v2.input_edges.append(self)
        following_edge.v1.output_edges.remove(following_edge)
        following_edge.v2.input_edges.remove(following_edge)
        self.seq += following_edge.seq[-1]
        self.cov = (self.cov * len(self) + following_edge.cov * len(following_edge)) / (len(self) + len(following_edge))

    def __str__(self):
        return self.seq if self.show_sequences else ''


class Vertex:
    show_sequences = False
    counter = itertools.count()

    def __init__(self, seq):
        self.seq = seq
        self.input_edges = []
        self.output_edges = []
        self.id = next(self.counter)

    def add_edge(self, other):
        for edge in other.input_edges:
            if self == edge.v1:
                edge.inc_coverage()
                break
        else:
            edge = Edge(self, other)
            other.input_edges.append(edge)
            self.output_edges.append(edge)

    def __str__(self):
        return self.seq if self.show_sequences else str(self.id)

    def can_compress(self):
        return len(self.input_edges) == len(self.output_edges) == 1


class Graph:
    k = None

    def __init__(self):
        self.all_vertices = {}

    def add_edge(self, seq1, seq2):
        if seq1 in self.all_vertices:
            v1 = self.all_vertices[seq1]
        else:
            v1 = self.all_vertices[seq1] = Vertex(seq1)

        if seq2 in self.all_vertices:
            v2 = self.all_vertices[seq2]
        else:
            v2 = self.all_vertices[seq2] = Vertex(seq2)

        v1.add_edge(v2)

    def add_seq(self, seq):
        for i in range(len(seq) - self.k):
            kmer_1 = seq[i:i + self.k]
            kmer_2 = seq[i + 1: i + self.k + 1]
            self.add_edge(kmer_1, kmer_2)

    def compress(self):
        to_delete = []
        for kmer, vertex in self.all_vertices.items():
            if vertex.can_compress():
                to_delete.append(kmer)
        for kmer in to_delete:
            self.all_vertices[kmer].input_edges[0].merge(self.all_vertices[kmer].output_edges[0])
            self.all_vertices.pop(kmer)

    def save_dot(self, outp):
        outp.write("digraph {\n")
        for v1 in self.all_vertices.values():
            for e in v1.output_edges:
                v2 = e.v2
                outp.write('{} -> {} [label="{}x{:.3f}"]\n'.format(
                    v1, v2, e, e.cov))
        outp.write("}")


complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def reverse_complement(seq):
    return ''.join(complement[nt] for nt in seq[::-1])


def read_fastq(f):
    for line in f:
        name = line.strip()
        seq = next(f).strip()
        next(f)
        next(f)
        yield name, seq


def read_fasta(f):
    name = None
    seq = None
    for line in f:
        if line.startswith('>'):
            if name:
                yield name, seq
            name = line.lstrip('>').strip()
            seq = ''
        else:
            seq += line.strip()
    yield name, seq


def read(f):
    if f.name.endswith('a'):
        return read_fasta(f)
    else:
        return read_fastq(f)


def main():
    parser = argparse.ArgumentParser(description='De Bruijn graph')
    parser.add_argument('-i', '--input', help='Input fastq', metavar='File',
                        type=argparse.FileType(), required=True)
    parser.add_argument('-k', help='k-mer size (default: 55)', metavar='Int',
                        type=int, default=55)
    parser.add_argument('-o', '--output', help='Output dot', metavar='File',
                        type=argparse.FileType('w'), required=True)
    parser.add_argument('-c', '--compress', help='Shrink graph', action='store_true')
    parser.add_argument('--vertex', help='Show vertex sequences', action='store_true')
    parser.add_argument('--edge', help='Show edge sequences', action='store_true')
    args = parser.parse_args()

    Graph.k = args.k
    Vertex.show_sequences = args.vertex
    Edge.show_sequences = args.edge

    graph = Graph()
    for name, seq in tqdm(read(args.input)):
        graph.add_seq(seq)
        graph.add_seq(reverse_complement(seq))

    print('Graph size: {}'.format(len(graph.all_vertices)))
    if args.compress:
        graph.compress()
        print('Compressed graph size: {}'.format(len(graph.all_vertices)))
    graph.save_dot(args.output)


if __name__ == '__main__':
    main()
