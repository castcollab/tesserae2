from collections import defaultdict

import attr
import contextlib
import io
import json
import sys

from pycortex.__main__ import main
import pycortex.test.builder as builder
import pycortex.test.runner as runner


@attr.s(slots=True)
class PycortexPrintOutputParser(object):
    output = attr.ib()

    def get_kmer_strings(self):
        return self.output.rstrip().split('\n')


@attr.s(slots=True)
class ViewExpectation(object):
    output = attr.ib()
    parser = attr.ib(init=False)
    kmer_strings = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.parser = PycortexPrintOutputParser(self.output)
        self.kmer_strings = self.parser.get_kmer_strings()

    def has_kmer(self, kmer_string):
        assert kmer_string in self.kmer_strings
        return self

    def has_n_kmers(self, n):
        assert len(self.kmer_strings) == n
        return self


class Test(object):
    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'AAG 1 ......G.',
            'ACC 1 .......T',
            'AGG 1 a......T',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()


class TestTermWithRecord(object):
    def test_prints_single_kmer(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('ACCAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CAA: CAA 1 1 .c...... ........'

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', 'CAA'])

        # then
        assert [expected_kmer] == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_one_missing_kmer(self, tmpdir):
        # given
        kmer_size = 3
        record = 'GGG'
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('AAAA')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmer = 'CCC: GGG 0 1 ........ ........'

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', record])

        # then
        assert [expected_kmer] == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers(self, tmpdir):
        # given
        record = 'ACCAA'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size).build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 1 ....A... ....A...',
            'CCA: CCA 1 1 a...A... a...A...',
            'CAA: CAA 1 1 .c...... .c......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_three_kmers_including_one_revcomp(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACC: ACC 1 1 .......T .......T',
            'AGG: CCT 1 1 A......t A......t',
            'AAG: CTT 1 1 .C...... .C......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', record])

        # then
        assert expected_kmers == PycortexPrintOutputParser(
            pycortex_output.getvalue()).get_kmer_strings()

    def test_prints_one_missing_and_one_revcomp_kmer(self, tmpdir):
        # given
        dna_sequence = 'ACCTT'
        search_record = 'ACTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(dna_sequence)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected_kmers = [
            'ACT: ACT 0 1 ........ .......T',
            'AAG: CTT 1 1 .C...... A.......',
        ]

        # when
        pycortex_output = io.StringIO()
        with contextlib.redirect_stdout(pycortex_output):
            main(['view', output_graph, '--record', search_record])
        expect = ViewExpectation(pycortex_output.getvalue())

        # then
        (expect
         .has_kmer(expected_kmers[0])
         .has_kmer(expected_kmers[1])
         .has_n_kmers(2))


def expect_zero_return_code(completed_process):
    stdout = completed_process.stdout.decode()
    if completed_process.returncode != 0:
        print(stdout)
        print(completed_process.stderr.decode(), file=sys.stderr)
        assert completed_process.returncode == 0


def listify_elements(iterable):
    return ([[e], e][int(isinstance(e, list))] for e in iterable)


@attr.s(slots=True)
class JsonNodeExpectation(object):
    node = attr.ib()
    n_colors = attr.ib(None)

    def has_coverages(self, *coverages):
        coverages = list(listify_elements(coverages))
        for coverage in coverages:
            assert len(coverage) == self.n_colors
        assert self.node['coverage'] == coverages
        return self


@attr.s(slots=True)
class JsonGraphExpectation(object):
    graph = attr.ib()
    colors = attr.ib(init=False)
    n_colors = attr.ib(init=False)
    node_id_by_repr = attr.ib(init=False)

    def __attrs_post_init__(self):
        print('With JSON graph: {}'.format(self.graph))
        assert sum(['is_missing' in n for n in self.graph['nodes']]) == 0
        assert sum(['is_missing' in e for e in self.graph['edges']]) == 0
        self.colors = self.graph['graph']['colors']
        self.n_colors = len(self.colors)
        self.node_id_by_repr = defaultdict(list)
        for node_id, node in enumerate(self.graph['nodes']):
            self.node_id_by_repr[node['repr']].append(node_id)

    def is_directed(self):
        assert self.graph['directed']
        return self

    def has_colors(self, colors):
        assert self.colors == colors
        return self

    def has_n_nodes(self, n):
        assert len(self.graph['nodes']) == n
        return self

    def has_node_repr(self, repr):
        nodes = list(filter(lambda n: n['repr'] == repr, self.graph['nodes']))
        assert len(nodes) > 0
        return JsonNodeExpectation(nodes[0], self.n_colors)

    def has_n_edges(self, n):
        assert len(self.graph['edges']) == n
        return self

    def has_repr_edge(self, source_repr, target_repr, color):
        source_id_list = self.node_id_by_repr[source_repr]
        target_id_list = self.node_id_by_repr[target_repr]
        assert len(source_id_list) == 1
        assert len(target_id_list) == 1
        print('source_repr={}, target_repr={}'.format(source_repr, target_repr))
        return self.has_edge(source_id_list[0], target_id_list[0], color)

    def has_edge(self, source, target, color):
        num_matching_edges = 0
        matching_edge = None
        for e in self.graph['edges']:
            edge = (e['source'], e['target'], e['key'])
            if edge == (source, target, color):
                num_matching_edges += 1
                matching_edge = edge
        assert matching_edge == (source, target, color)
        assert num_matching_edges == 1
        return self


class TestOutputTypeJSON(object):
    def test_outputs_json(self, tmpdir):
        # given
        record = 'ACCTT'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        # when
        completed_process = (runner
                             .Pycortex(True)
                             .view(['--record', record, '--output-type', 'json', output_graph]))
        stdout = completed_process.stdout.decode()

        # then
        assert completed_process.returncode == 0, completed_process
        expect = JsonGraphExpectation(json.loads(stdout))
        expect.is_directed()

    def test_collapse_kmer_unitigs_option(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record1)
                        .with_dna_sequence(record2)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = (runner
                             .Pycortex(True)
                             .view(['--record', record1,
                                    '--output-type', 'json',
                                    output_graph]))

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout.decode()
        expect = JsonGraphExpectation(json.loads(stdout))

        expect.has_colors([0, 1])
        expect.has_n_nodes(3)
        expect.has_node_repr('AAACC').has_coverages([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('C').has_coverages([1, 1])
        expect.has_node_repr('GAA').has_coverages([2, 1], [1, 1], [1, 1])

        for edge in [(0, 1, 0), (0, 1, 1), (1, 2, 0), (1, 2, 1), (0, 2, 0)]:
            expect.has_edge(*edge)
        expect.has_n_edges(5)

    def test_collapse_kmer_unitigs_option_with_missing_kmers(self, tmpdir):
        # given
        record1 = 'AAACCCGAA'
        record2 = 'ACCG'
        query_record = record1 + 'G'
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence(record1)
                        .with_dna_sequence(record2)
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))
        runner.Mccortex(kmer_size).view(output_graph)

        # when
        completed_process = (runner
                             .Pycortex(True)
                             .view(['--record', query_record,
                                    '--output-type', 'json',
                                    output_graph]))

        # then
        expect_zero_return_code(completed_process)

        stdout = completed_process.stdout.decode()
        expect = JsonGraphExpectation(json.loads(stdout))

        expect.has_n_nodes(4)
        expect.has_node_repr('AAACC').has_coverages([1, 1], [1, 1], [2, 1])
        expect.has_node_repr('C').has_coverages([1, 1])
        expect.has_node_repr('GAA').has_coverages([2, 1], [1, 1], [1, 1])
        expect.has_node_repr('G').has_coverages([0, 1])

        for color in [0, 1]:
            for edge in [['AAACC', 'C'], ['C', 'GAA']]:
                expect.has_repr_edge(edge[0], edge[1], color)
        expect.has_repr_edge('GAA', 'G', 1)
        expect.has_repr_edge('AAACC', 'GAA', 0)
        expect.has_n_edges(6)
