import pandas as pd


class BedFile(object):
    """Class definition for a BED file (contains a list of intervals)."""

    def __init__(self, x):
        """Initialize BedFile with string or pandas DataFrame."""
        if type(x) is str:
            with open(x) as f:
                rows = (line.split() for line in f)
                self.intervals = [Interval(row[0], int(row[1]), int(row[2]))
                                  for row in rows]
        elif type(x) is pd.DataFrame:
            # assumes column chr/start/end names
            self.intervals = [Interval(row['chr'], row['start'], row['end'])
                              for i, row in x.iterrows()]


class Interval(object):
    """Class definition for a BED interval."""

    def __init__(self, chrom, start, end, **kwargs):
        """
        Attributes
        ----------
        chrom : str
        start : int
        end : int
        """
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = end - start
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def overlaps(self, interval):
        """Indicates whether there is an overlap with another Interval."""
        return ((self.chrom == interval.chrom) and
                (self.start < interval.end) and
                (interval.start < self.end))

    def __hash__(self):
        """Hash as (chrom, start, end)"""
        return hash((self.chrom, self.start, self.end))

    def __repr__(self):
        """Represent as 'Interval(chrom, start, end)'"""
        return 'Interval({}, {}, {})'.format(self.chrom, self.start, self.end)
