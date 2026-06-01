import csv

vcf_in = snakemake.input[0]
csv_out = snakemake.output.csv
txt_out = snakemake.output.txt

csq_fields = []
with open(vcf_in) as f:
    for line in f:
        if line.startswith("##INFO=<ID=CSQ"):
            fmt_str = line.split("Format: ")[1].strip().rstrip('">').strip()
            csq_fields = fmt_str.split("|")
            break

base_header = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","WES","WGS"]
full_header = base_header + csq_fields

with open(vcf_in) as vcf, \
     open(csv_out, "w", newline="") as csvf, \
     open(txt_out, "w") as txtf:

    writer = csv.writer(csvf)
    writer.writerow(full_header)
    txtf.write("\t".join(full_header) + "\n")

    for line in vcf:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue

        chrom,pos,id_,ref,alt,qual,filter_,info,fmt = fields[:9]
        wes, wgs = fields[9], fields[10]

        if filter_ != "PASS": continue
        if "missense_variant" not in info: continue

        # Parse CSQ first before filtering on SIFT
        csq_values = [""] * len(csq_fields)
        for info_field in info.split(";"):
            if info_field.startswith("CSQ="):
                first_transcript = info_field[4:].split(",")[0]
                csq_values = first_transcript.split("|")
                break

        # Get SIFT from parsed CSQ fields
        sift_idx = csq_fields.index("SIFT")
        sift_val = csq_values[sift_idx] if len(csq_values) > sift_idx else ""

        # Strict SIFT filter - excludes empty, tolerated, deleterious_low_confidence
        if not sift_val.startswith("deleterious(0)"): continue

        # Genotype filter
        wes_gt = wes.split(":")[0]
        wgs_gt = wgs.split(":")[0]
        if wes_gt != "1/1" or wgs_gt != "1/1": continue

        # Alt read depth filter
        wes_ad = wes.split(":")[1].split(",")
        wgs_ad = wgs.split(":")[1].split(",")
        if int(wes_ad[1]) < 10 or int(wgs_ad[1]) < 10: continue

        row = [chrom,pos,id_,ref,alt,qual,filter_,fmt,wes,wgs] + csq_values
        writer.writerow(row)
        txtf.write("\t".join(row) + "\n")