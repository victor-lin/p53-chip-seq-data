import re


class SeqSample(object):
    def __init__(self, sample_name):
        """
        Attributes
        ----------
        sample_name : str
        cell_type : str
        treatment_type : str
        treatment_time : int
        treatment_repeat : str
        """
        self.sample_name = sample_name
        self.__parse_sample_name()

    def __parse_sample_name(self):
        """Parse the sample name and set attributes."""
        pattern = '(.*)(P53)(XR|NT)(\d+)([A-Z]?|Ctr)?.*'
        vals = re.findall(pattern, self.sample_name.replace('_', ''))[0]
        self.cell_type = vals[0]
        self.treatment_type = vals[2]
        self.treatment_time = vals[3]
        if vals[3]:
            self.treatment_repeat = vals[4]

    def __str__(self):
        slots = ['cell_type', 'treatment_type', 'treatment_time']
        slot_val_strs = ['{}={}'.format(k, getattr(self, k))
                         for k in slots if hasattr(self, k)]
        return 'SeqSample({})'.format(','.join(slot_val_strs))
