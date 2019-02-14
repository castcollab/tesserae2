from cortexpy.tesserae import Tesserae


class TestTesserae:
    def test_mosaic_alignment_on_short_query_and_two_templates(self):
        # given
        query = "GTAGGCGAGATGACGCCAT"
        targets = ["GTAGGCGAGTCCCGTTTATA", "CCACAGAAGATGACGCCATT"]

        # when
        t = Tesserae()
        p = t.align(query, targets)

        # then
        assert len(p) == 3

        assert p[1][0] == 'template0'
        assert p[1][1] == 'GTAGGCG'
        assert p[1][2] == 0
        assert p[1][3] == 6

        assert p[2][0] == 'template1'
        assert p[2][1] == '       AGATGACGCCAT'
        assert p[2][2] == 7
        assert p[2][3] == 18
