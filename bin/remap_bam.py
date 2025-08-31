#!/usr/bin/env python3
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

        header_dict["SQ"] = [{"SN": rname, "LN": ln} for rname in rnames]   # new references

        with pysam.AlignmentFile("temp.bam", "wb", header=header_dict) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                qname = read.query_name
                if qname in qname_to_rname:
                    bam_out.write(read)

    # second pass
    with pysam.AlignmentFile("temp.bam", "rb") as bam_temp:
        with pysam.AlignmentFile(output_bam, "wb", header=bam_temp.header) as bam_out:
            for read in bam_temp.fetch(until_eof=True):
                qname = read.query_name
                if qname in qname_to_rname:
                    rname = qname_to_rname[qname]
                    rid = bam_out.get_tid(rname)
                    read.reference_id = rid
                    bam_out.write(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--barcode_map", required=True, type=str)
    parser.add_argument("-i", "--input_bam", required=True, type=str)
    parser.add_argument("-o", "--output_bam", required=True, type=str)
    args = parser.parse_args()

    main(args.barcode_map, args.input_bam, args.output_bam)
