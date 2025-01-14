from io import StringIO, BytesIO
import tarfile
import zipfile

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from libraries.multilift.liftover import annotate, Lifter, liftover
from libraries.multilift.msa import align, generate_consensus
from libraries.multilift.utils import add_to_archive, sniff_filetype


def run_multilift(genomes, sequences, aligner, output_format):
    """
    Run multilift for a set of genomes and sequences.

    Args:
        genomes (list): List of genome names.
        sequences (dict): Dictionary mapping (genome, file, seq_id) to SeqRecord.
        aligner (str): Aligner to use for sequence alignment.
        output_format (str): Output archive format (e.g., ".tar.gz" or ".zip").

    Returns:
        BytesIO: Archive containing multilift results.
    """
    multilift_download = BytesIO()
    maf_alignments = []
    igv_genomes = []
    L = Lifter()

    with (
        tarfile.open(fileobj=multilift_download, mode="w:gz")
        if output_format == ".tar.gz"
        else zipfile.ZipFile(multilift_download, "w") as Arc
    ):
        # Perform alignment for each genome
        for genome in genomes:
            group_sequences = [
                seq for (g, _, _), seq in sequences.items() if g == genome
            ]
            if group_sequences:
                with StringIO() as F:
                    SeqIO.write(group_sequences, F, "fasta")
                    returncode, result = align(F, aligner)
                if returncode:
                    raise RuntimeError(
                        f"Error making {aligner} alignment for genome: {genome}"
                    )

                aln = AlignIO.read(result, "fasta")
                L.add_alignment(aln, genome)

                # Write alignment as fasta
                add_to_archive(
                    Arc,
                    BytesIO(bytes(result.getvalue(), "utf-8")),
                    f"{genome}/alignment.fa",
                    output_format,
                )

                # Store consensus
                cons = Seq(generate_consensus(aln))
                igv_genomes.append(SeqRecord(id=genome, description="", seq=cons))

                # Add to MAF alignments
                maf_alignments.append(
                    MultipleSeqAlignment(
                        [
                            SeqRecord(
                                id=f"multilift.{genome}", description="", seq=cons
                            ),
                            *aln,
                        ]
                    )
                )

        # Write consensus sequences
        with StringIO() as F:
            SeqIO.write(igv_genomes, F, "fasta")
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), "utf-8")),
                "genome/multilift_consensus.fa",
                output_format,
            )

        # Write MAF alignment
        with StringIO() as F:
            AlignIO.write(maf_alignments, F, "maf")
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), "utf-8")),
                "alignment/multilift.maf",
                output_format,
            )

    multilift_download.seek(0)
    return multilift_download


def perform_liftover(genomes, sequences, consensus, output_format):
    """
    Perform liftover for a set of genomes and sequences.

    Args:
        genomes (list): List of genome names.
        sequences (dict): Dictionary mapping (genome, file, seq_id) to SeqRecord.
        consensus (bool): Whether to use consensus sequences for liftover.
        output_format (str): Output archive format (e.g., ".tar.gz" or ".zip").

    Returns:
        BytesIO: Archive containing liftover results.
    """
    liftover_download = BytesIO()
    L = Lifter()

    with (
        tarfile.open(fileobj=liftover_download, mode="w:gz")
        if output_format == ".tar.gz"
        else zipfile.ZipFile(liftover_download, "w") as Arc
    ):
        for genome in genomes:
            genome_sequences = {k: v for k, v in sequences.items() if k[0] == genome}

            for (genome_name, file_name, seq_id), seq in genome_sequences.items():
                if consensus:
                    new_ext, lift_file = annotate(
                        StringIO(str(seq.seq)), L, genome_sequences, genome_name
                    )
                else:
                    new_ext, lift_file = liftover(
                        StringIO(str(seq.seq)), "fasta", L, genome_name
                    )

                add_to_archive(
                    Arc,
                    BytesIO(bytes(lift_file.getvalue(), "utf-8")),
                    f"liftover/{genome}/{file_name}{new_ext}",
                    output_format,
                )

    liftover_download.seek(0)
    return liftover_download


def parse_sequence(genome, sequence, multilift_sequences):
    """
    Parse sequences from a file or string.

    Args:
        genome (str): Genome name.
        sequence (str): Sequence data.
        multilift_sequences (dict): Dictionary mapping (genome, file, seq_id) to SeqRecord

    Returns:
        dict: Dictionary mapping (genome, file, seq_id, group) to SeqRecord.
    """
    ftype, _ = sniff_filetype(sequence.name)
    try:
        with StringIO(sequence.read().decode("utf-8")) as F:
            multilift_sequences.update(
                {(genome, sequence.name, seq.id): seq for seq in SeqIO.parse(F, ftype)}
            )
            return True
    except Exception as e:
        print("Failed to parse sequence file:", e)
        return False
