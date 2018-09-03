from scripts.binding_dataset.bed import Interval


def test_bed_module():
    assert Interval('chr2L', 5, 11).overlaps(Interval('chr2L', 10, 20))
    assert not Interval('chr2L', 5, 10).overlaps(Interval('chr2L', 10, 20))
    assert Interval('chr2L', 19, 30).overlaps(Interval('chr2L', 10, 20))
    assert not Interval('chr2L', 20, 30).overlaps(Interval('chr2L', 10, 20))
    assert Interval('chr2L', 11, 19).overlaps(Interval('chr2L', 10, 20))
    assert Interval('chr2L', 5, 30).overlaps(Interval('chr2L', 10, 20))
