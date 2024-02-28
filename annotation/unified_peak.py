

"""
This script defines classes and functions for working with peaks and chromosomes.
It provides functionality to read peak and chromosome information from files,
assign peaks to chromosomes based on their chromosomal location, and remove regions without peaks.
"""


class Peak:
    def __init__(self, chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak):
        """
        Initialize a Peak object.

        Args:
            chrom (str): The chromosome of the peak.
            start (int): The start position of the peak.
            end (int): The end position of the peak.
            name (str): The name of the peak.
            score (int): The score of the peak.
            strand (str): The strand of the peak.
            signalValue (float): The signal value of the peak.
            pValue (float): The p-value of the peak.
            qValue (float): The q-value of the peak.
            peak (int): The peak position relative to the start position.

        Returns:
            None
        """
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.signalValue = float(signalValue)
        self.pValue = float(pValue)
        self.qValue = float(qValue)
        self.peak = int(peak)

    def peak_position(self):
        """
        Calculate the absolute position of the peak.

        Returns:
            int: The absolute position of the peak.
        """
        return self.start + self.peak
    
    def __str__(self):
        """
        Return a string representation of the Peak object.

        Returns:
            str: A string representation of the Peak object.
        """
        return f"Peak {self.name} at {self.chrom}:{self.start}-{self.end} with score {self.score} and peak position {self.peak_position()}"
    
    __repr__ = __str__

    @staticmethod
    def read_peaks(file_path):
        """
        Read peaks from a file and return a list of Peak objects.

        Args:
            file_path (str): The path to the file containing the peaks.

        Returns:
            list: A list of Peak objects.
        """
        peaks = []
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                peaks.append(Peak(*parts))
        return peaks


class Chromosome:
    def __init__(self, name: str, length: int, region_size: int = 100):
        """
        Initialize a Chromosome object.

        Args:
            name (str): The name of the chromosome.
            length (int): The length of the chromosome.
            region_size (int, optional): The size of each region. Defaults to 100.
        """
        self.name = name
        self.region_size = region_size
        self.intervals = {i: [] for i in range(0, length, region_size)}

    def add_peak(self, peak):
        """
        Add a peak to the chromosome.

        Args:
            peak (Peak): The peak object to be added.
        """
        interval_index = peak.peak_position() // self.region_size * self.region_size
        if interval_index in self.intervals:
            self.intervals[interval_index].append(peak)
            #print(f"Added peak {peak.name} to interval {interval_index} in chromosome {self.name}")
        else:
            print(f"Peak {peak.name} is out of range for chromosome {self.name}")

    def __str__(self):
        """
        Return a string representation of the chromosome.

        Returns:
            str: The string representation of the chromosome.
        """
        return f"Chromosome {self.name} with {len(self.intervals)} intervals "

    __repr__ = __str__

    @staticmethod
    def read_chromosomes(file_path: str, region_size: int = 100):
        """
        Read chromosomes from a file.

        Args:
            file_path (str): The path to the file containing chromosome information.
            region_size (int, optional): The size of each region. Defaults to 100.

        Returns:
            dict: A dictionary of Chromosome objects, where the keys are chromosome names.
        """
        chromosomes = {}
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                chrom_name, length = parts[0], int(parts[1])
                chromosomes[chrom_name] = Chromosome(chrom_name, length, region_size)
        return chromosomes


def assign_peaks_to_chromosomes(peaks, chromosomes):
    """
    Assigns peaks to chromosomes based on their chromosomal location.

    Args:
        peaks (list): List of peak objects.
        chromosomes (dict): Dictionary containing chromosome objects.

    Returns:
        None
    """
    for peak in peaks:
        if peak.chrom in chromosomes:
            chromosomes[peak.chrom].add_peak(peak)
        else:
            print(f"Chromosome {peak.chrom} not found.")

def remove_region_without_peak(chromosomes):
    """
    Remove regions without peaks from the given chromosomes.

    Args:
        chromosomes (dict): A dictionary containing chromosomes and their intervals.

    Returns:
        None
    """
    for chrom in chromosomes.values():
        intervels_to_remove = []
        for interval in list(chrom.intervals.keys()):
            # first record the interval
            if len(chrom.intervals[interval]) == 0:
                intervels_to_remove.append(interval)
        # then remove the interval
        for interval in intervels_to_remove:
            del chrom.intervals[interval]

        
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assign peaks to chromosomes.")
    parser.add_argument("-r", "--region_size", help="The size of each region.", default=100, type=int)
    parser.add_argument("peak_dir", help="Path to the file containing peaks.")
    parser.add_argument("chromosome_file", help="Path to the file containing chromosome information.")
    parser.add_argument("out_peak_bed", help="Path to the output file.")
    args = parser.parse_args()


    # get all peak in /data5/xzg_data/GRN/Plant_Unified_Peak/peak which end with .narrowPeak

    import os
    peaks = []
    for root, dirs, files in os.walk(args.peak_dir):
        for file in files:
            if file.endswith(".narrowPeak"):
                peaks.extend(Peak.read_peaks(os.path.join(root, file)))

    chromosome = args.chromosome_file
    # get all chromosome in /data5/xzg_data/GRN/Plant_Unified_Peak/chromosomes
    chromosomes = Chromosome.read_chromosomes(chromosome, args.region_size)

    assign_peaks_to_chromosomes(peaks, chromosomes)
    remove_region_without_peak(chromosomes)
    # save the result in bed format in /data5/xzg_data/GRN/Plant_Unified_Peak/peak_bed
    with open(args.out_peak_bed, "w") as f:
        for chrom in chromosomes.values():
            for interval in chrom.intervals:
                # if the peak num < 2 continue
                #if len(chrom.intervals[interval]) < 2:
                #    continue
                # save the interval: chrom, start, end, peak_name, peak_num strand(.)
                f.write(f"{chrom.name}\t{interval}\t{interval + chrom.region_size}\t{chrom.name}_{interval}\t{len(chrom.intervals[interval])}\t.\n")
