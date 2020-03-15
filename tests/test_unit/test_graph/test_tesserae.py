from cortexpy.tesserae import Tesserae
from numpy import sqrt


class TestTesserae:
    def __init__(self):
        # given
        self.query = "GTAGGCGAGATGACGCCAT"
        self.targets = ["GTAGGCGAGTCCCGTTTATA", "CCACAGAAGATGACGCCATT"]
        # when
        self.t = Tesserae(mem_limit=False)
        self.p = self.t.align(self.query, self.targets)

    def test_mosaic_alignment_on_short_query_and_two_templates(self):
        # then
        assert len(self.p) == 3

        assert self.p[1][0] == 'template0'
        assert self.p[1][1] == 'GTAGGCG'
        assert self.p[1][2] == 0
        assert self.p[1][3] == 6

        assert self.p[2][0] == 'template1'
        assert self.p[2][1] == '       AGATGACGCCAT'
        assert self.p[2][2] == 7
        assert self.p[2][3] == 18

        assert self.t.llk >= -27
        assert self.t.llk <= -26

    def test_mosaic_alignment_on_short_query_and_two_templates_reduced_memory_usage(self):
        # when
        t = Tesserae(mem_limit=True)
        p = t.align(self.query, self.targets)
        max_mem_limit = sqrt(len(self.query)) + 1

        # then
        assert p == self.p

        assert t.llk >= -27
        assert t.llk <= -26

        assert t.traceback_limit <= max_mem_limit
        assert t.states_to_save <= max_mem_limit
