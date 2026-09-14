import gzip
import logging
import os


INPUT_FILE = "Data/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38.txt.gz"
OUTPUT_FILE = "Data/databases/haplotypes/All_indep_gwas_sumstats_AI_hg38_allele_matches.txt.gz"
PROGRESS_INTERVAL = 1_000_000

CHROM_INDEX = 2
POSITION_INDEX = 3
EA_INDEX = 5
NEA_INDEX = 6
EAF_INDEX = 7
BETA_INDEX = 8

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_eaf(value):
    """Return a valid allele frequency as a float, otherwise None."""
    try:
        eaf = float(value)
    except (TypeError, ValueError):
        return None
    return eaf if 0 <= eaf <= 1 else None


def unordered_pair(ea, nea):
    """Treat EA/NEA and NEA/EA as the same allele pair."""
    return tuple(sorted((ea, nea)))


def identify_retained_alleles(input_file):
    """First pass: choose the majority allele pair and ordering at each locus."""
    # Most loci contain one unordered allele pair. Store that common case
    # compactly and allocate an alternatives dictionary only for loci where a
    # genuinely different pair occurs.
    locus_states = {}
    alternative_states = {}
    total_records = 0

    logger.info("Pass 1/2: identifying the majority allele pair at each locus")
    with gzip.open(input_file, "rt") as source:
        next(source)
        for raw_line in source:
            fields = raw_line.rstrip("\n").split("\t")
            total_records += 1
            locus = (fields[CHROM_INDEX], fields[POSITION_INDEX])
            orientation = (fields[EA_INDEX], fields[NEA_INDEX])
            pair = unordered_pair(*orientation)
            eaf = parse_eaf(fields[EAF_INDEX])

            if locus not in locus_states:
                # pair, count, first ordered representation, first usable EAF
                locus_states[locus] = [pair, 1, orientation, eaf]
            elif pair == locus_states[locus][0]:
                locus_states[locus][1] += 1
                if locus_states[locus][3] is None and orientation == locus_states[locus][2]:
                    locus_states[locus][3] = eaf
            else:
                alternatives = alternative_states.setdefault(locus, {})
                if pair not in alternatives:
                    alternatives[pair] = [1, orientation, eaf]
                else:
                    alternatives[pair][0] += 1
                    if alternatives[pair][2] is None and orientation == alternatives[pair][1]:
                        alternatives[pair][2] = eaf

            if total_records % PROGRESS_INTERVAL == 0:
                logger.info("Pass 1/2: processed %s records", f"{total_records:,}")

    tied_pair_loci = 0
    expected_removed_records = 0
    logger.info(
        "Pass 1/2: resolving majority votes at %s loci with multiple allele pairs",
        f"{len(alternative_states):,}",
    )
    for locus, state in locus_states.items():
        candidates = [(state[0], state[1], state[2], state[3])]
        for pair, alternative in alternative_states.get(locus, {}).items():
            candidates.append((pair, alternative[0], alternative[1], alternative[2]))

        highest_count = max(candidate[1] for candidate in candidates)
        winners = [candidate for candidate in candidates if candidate[1] == highest_count]
        winner = winners[0]
        if len(winners) > 1:
            tied_pair_loci += 1

        expected_removed_records += sum(candidate[1] for candidate in candidates) - winner[1]
        # Retain only the chosen ordered representation and its reference EAF.
        locus_states[locus] = (*winner[2], winner[3])

    loci_with_multiple_pairs = len(alternative_states)
    logger.info(
        "Pass 1/2 complete: %s records across %s loci; %s loci have multiple allele pairs",
        f"{total_records:,}",
        f"{len(locus_states):,}",
        f"{loci_with_multiple_pairs:,}",
    )
    if tied_pair_loci:
        logger.warning(
            "%s loci have a tie between majority allele pairs; the first encountered tied pair will be retained",
            f"{tied_pair_loci:,}",
        )

    return locus_states, {
        "input_records": total_records,
        "unique_loci": len(locus_states),
        "loci_with_multiple_pairs": loci_with_multiple_pairs,
        "tied_pair_loci": tied_pair_loci,
        "expected_removed_records": expected_removed_records,
    }


def write_harmonized_file(input_file, output_file, retained, pass_one_stats):
    """Second pass: retain, flip, or remove each association record."""
    temporary_output = f"{output_file}.tmp"
    unchanged_records = 0
    flipped_records = 0
    removed_records = 0
    flipped_loci = set()
    removed_loci = set()
    eaf_opposite = 0
    eaf_not_opposite = 0
    eaf_missing = 0

    logger.info("Pass 2/2: writing harmonized associations to %s", output_file)
    try:
        with gzip.open(input_file, "rt") as source, gzip.open(temporary_output, "wt") as destination:
            destination.write(next(source))

            for record_number, raw_line in enumerate(source, start=1):
                fields = raw_line.rstrip("\n").split("\t")
                locus = (fields[CHROM_INDEX], fields[POSITION_INDEX])
                retained_ea, retained_nea, reference_eaf = retained[locus]
                retained_pair = unordered_pair(retained_ea, retained_nea)
                record_pair = unordered_pair(fields[EA_INDEX], fields[NEA_INDEX])

                if record_pair != retained_pair:
                    removed_records += 1
                    removed_loci.add(locus)
                    continue

                if (fields[EA_INDEX], fields[NEA_INDEX]) == (retained_ea, retained_nea):
                    unchanged_records += 1
                else:
                    record_eaf = parse_eaf(fields[EAF_INDEX])
                    if reference_eaf is None or record_eaf is None:
                        eaf_missing += 1
                    elif ((reference_eaf < 0.5 and record_eaf > 0.5) or
                          (reference_eaf > 0.5 and record_eaf < 0.5)):
                        eaf_opposite += 1
                    else:
                        eaf_not_opposite += 1

                    fields[EA_INDEX], fields[NEA_INDEX] = retained_ea, retained_nea
                    fields[BETA_INDEX] = str(-float(fields[BETA_INDEX]))
                    if record_eaf is not None:
                        fields[EAF_INDEX] = str(1.0 - record_eaf)
                    flipped_records += 1
                    flipped_loci.add(locus)

                destination.write("\t".join(fields) + "\n")

                if record_number % PROGRESS_INTERVAL == 0:
                    logger.info("Pass 2/2: processed %s records", f"{record_number:,}")

        os.replace(temporary_output, output_file)
    except Exception:
        if os.path.exists(temporary_output):
            os.remove(temporary_output)
        raise

    output_records = unchanged_records + flipped_records
    if removed_records != pass_one_stats["expected_removed_records"]:
        raise RuntimeError(
            "Removal count differs between passes: "
            f"expected {pass_one_stats['expected_removed_records']}, observed {removed_records}"
        )
    if output_records + removed_records != pass_one_stats["input_records"]:
        raise RuntimeError("Output and removed record counts do not reconcile with the input")

    logger.info("Harmonization complete")
    logger.info("Input association records: %s", f"{pass_one_stats['input_records']:,}")
    logger.info("Output association records: %s", f"{output_records:,}")
    logger.info("Unchanged records: %s", f"{unchanged_records:,}")
    logger.info("Flipped records: %s across %s loci", f"{flipped_records:,}", f"{len(flipped_loci):,}")
    logger.info("  Flips with EAFs on opposite sides of 0.5: %s", f"{eaf_opposite:,}")
    logger.info("  Flips with EAFs not on opposite sides of 0.5: %s", f"{eaf_not_opposite:,}")
    logger.info("  Flips with unavailable EAF: %s", f"{eaf_missing:,}")
    logger.info("Removed minority-pair records: %s across %s loci", f"{removed_records:,}", f"{len(removed_loci):,}")
    logger.info("Majority-pair ties: %s loci", f"{pass_one_stats['tied_pair_loci']:,}")


def main():
    retained, pass_one_stats = identify_retained_alleles(INPUT_FILE)
    write_harmonized_file(INPUT_FILE, OUTPUT_FILE, retained, pass_one_stats)


if __name__ == "__main__":
    main()
