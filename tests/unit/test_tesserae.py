from tesserae import Tesserae
from numpy import sqrt
import pytest


class TesseraeFixture:
    def __init__(self):
        self.query = "GTAGGCGAGATGACGCCAT"
        self.targets = ["GTAGGCGAGTCCCGTTTATA", "CCACAGAAGATGACGCCATT"]
        # when
        self.t = Tesserae(mem_limit=False)
        self.p = self.t.align(self.query, self.targets)


@pytest.fixture
def fixture():
    return TesseraeFixture()


class TestTesserae:
    def test_mosaic_alignment_on_short_query_and_two_templates(self, fixture):

        # then
        assert len(fixture.p) == 3

        assert fixture.p[1][0] == "template0"
        assert fixture.p[1][1] == "GTAGGCG"
        assert fixture.p[1][2] == 0
        assert fixture.p[1][3] == 6

        assert fixture.p[2][0] == "template1"
        assert fixture.p[2][1] == "       AGATGACGCCAT"
        assert fixture.p[2][2] == 7
        assert fixture.p[2][3] == 18

        assert fixture.t.llk >= -27
        assert fixture.t.llk <= -26

    def test_mosaic_alignment_on_short_query_and_two_templates_reduced_memory_usage(
        self, fixture
    ):
        # when
        t = Tesserae(mem_limit=True)
        p = t.align(fixture.query, fixture.targets)
        max_mem_limit = sqrt(len(fixture.query)) + 1

        # then
        assert p == fixture.p

        assert t.llk >= -27
        assert t.llk <= -26

        assert t.traceback_limit <= max_mem_limit
        assert t.states_to_save <= max_mem_limit
