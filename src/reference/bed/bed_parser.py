"""Based on 'Enhancing GNN training for de novo genome assembly by targeting decision nodes'."""

import sys
from collections import defaultdict
from pathlib import Path


def parse_bed_file(bed_file_path: Path) -> dict[str, list]:
    """Parse the bed file to extract regions of interest.

    Bed file contains coluumn
    Args:
        bed_file_path: Path to file containing bed entries, notably centromere locations
    Return:
        dict containing mapping between chromosome name and all detected regions (start, stop)
    """
    bed_entries = defaultdict(list)
    with open(bed_file_path) as bed:
        for line in bed:
            columns = line.strip().split()
            if columns[0] in ["#", "track", "browser"]:
                # sometimes there is a header in file so skip this lines first
                continue
            if len(columns) >= 3:
                chromosome = columns[0]
                try:
                    bed_start = int(columns[1])
                    # the stop is 1-index but we don't care about that
                    bed_stop = int(columns[2])
                except ValueError as e:
                    print(f"{e} error occurred, skipping line {line}")
                    continue
                bed_entries[chromosome].append((bed_start, bed_stop))

    return bed_entries


def merge_intervals(intervals: list[tuple], gap: int) -> list[tuple]:
    """Either concatenates overlapping ranges (or one that are less then `gap` nucleotides apart) or simply appends
    intervals to list if they are to far apart."""
    # Sort intervals based on the start position
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        # If the list of merged intervals is empty or if the current interval does not overlap
        # with the previous one (considering a gap), simply append it.
        if not merged or merged[-1][1] + gap < interval[0]:
            merged.append(interval)
        else:
            # There is an overlap or a small gap, so extend the previous interval with the current one.
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))

    return merged


def merge_bed_regions(bed_intervals: dict, neighborhood: int = 100000, gap: int = 20000) -> dict[str, tuple]:
    """Takes a per chromosome regions from bed and tries to merge them into a single contiguous region.

    This is achieved by connecting consecutive regions which are no
    more then `gap` nucleotides apart. The largest region after concatenating is additionally
    expanded left and right by `neighborhood` nucleotides
    """
    largest_regions = {}  # Dictionary to store the result
    for current_chr in bed_intervals.keys():
        # Merge the intervals for the current chromosome
        merged_intervals = merge_intervals(bed_intervals[current_chr], gap=gap)
        # Find the largest interval
        if merged_intervals:  # Ensure there is at least one interval to avoid errors
            largest_interval = max(merged_intervals, key=lambda x: x[1] - x[0])
            # Adjust the interval by subtracting and adding max_distance to the start and end, respectively
            start_interval = largest_interval[0] - neighborhood
            end_interval = largest_interval[1] + neighborhood
            # Store the adjusted interval in the dictionary
            largest_regions[current_chr] = (max(0, start_interval), end_interval)
        else:
            # In case there are no centromeric intervals for the chromosome, you could choose to either skip it
            # or assign a default value such as None or an empty tuple.
            largest_regions[current_chr] = None  # or use (0, 0) or similar based on your requirements

    print(f"Largest regions {len(largest_regions)}")
    return largest_regions


def pretty_print_regions(regions: dict, output_file: Path):
    """Store the selected regions in chromosome to a file.

    Outputs chromosome name, length of the range start and end offsets in the reference chromosome
    """
    with open(output_file, "w") as f:
        header = f"{'Name'.ljust(15)}\t{'Lenght'.rjust(10)}\t{'Start'.rjust(10)}\t{'Stop'.rjust(10)}\n"
        f.write(header)
        for chr, subrange in regions.items():
            length = str(subrange[1] - subrange[0])
            start, stop = str(subrange[0]), str(subrange[1])
            line = f"{chr.ljust(15)}\t{length.rjust(10)}\t{start.rjust(10)}\t{stop.rjust(10)}\n"
            f.write(line)


if __name__ == "__main__":
    bed_file_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    largest_regions = merge_bed_regions(parse_bed_file(bed_file_path))
    pretty_print_regions(largest_regions, output_path)
