#!/usr/bin/env python3
"""
Thank you, Lorenz (https://github.com/lorenzwidmer), for the idea with the BAM remapping!
"""
import argparse
import csv
import pysam


def main(barcode_map, input_bam, output_bam):
    with open(barcode_map, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)
        qname_to_rname = {row['read_id']: row['cluster_id'] for row in rows}  # dict mapping qname to rname
        rnames = {row['cluster_id'] for row in rows}                          # uniqe rnames

    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        header_dict = bam_in.header.to_dict()

        if len(header_dict["SQ"]) != 1:
            raise ValueError("Not exactly one reference found!")

        ln = header_dict["SQ"][0]["LN"]

        header_dict["SQ"] += [{"SN": rname, "LN": ln} for rname in rnames]   # new references
        header = pysam.AlignmentHeader.from_dict(header_dict)

    with pysam.AlignmentFile(input_bam, "rb", header=header) as bam_in:
        with pysam.AlignmentFile(output_bam, "wb", header=header) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                qname = read.query_name
                if qname in qname_to_rname:
                    rname = qname_to_rname[qname]

                    rid = header.get_tid(rname)
                    new_read = pysam.AlignedSegment.fromstring(read.to_string(), header)
                    new_read.reference_id = rid
                    
                    bam_out.write(new_read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--barcode_map", required=True, type=str)
    parser.add_argument("-i", "--input_bam", required=True, type=str)
    parser.add_argument("-o", "--output_bam", required=True, type=str)
    args = parser.parse_args()

    main(args.barcode_map, args.input_bam, args.output_bam)
