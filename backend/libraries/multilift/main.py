from io import StringIO, BytesIO
import tarfile
import zipfile

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from libraries.multilift.liftover import annotate, Lifter, liftover
from libraries.multilift.msa import align, generate_consensus
from libraries.multilift.utils import (
    add_to_archive,
    sniff_filetype,
)


def get_multilift_sequences(multilift_sequences, genome_file, genome):
    ftype, _ = sniff_filetype(genome_file.name)
    try:
        with StringIO(genome_file.read().decode("utf-8")) as F:
            multilift_sequences.update(
                {
                    (genome, genome_file.name, seq.id, None): seq
                    for seq in SeqIO.parse(F, ftype)
                }
            )
    except:
        print(f"Error reading sequences: {genome_file.name}", 3)

    return multilift_sequences


def parse_sequence_format(sequences: dict):
    """
    Parses sequences from the format provided to a dictionary

    Args:
    - sequences: Dictionary from key to group name

    Returns:
    - parsed_sequences: Dictionary mapping from tuple of (genome, fname, seqid) to group

    """
    res = {}
    for key in sequences.keys():
        new_key = tuple(key.split(","))
        res[new_key] = sequences[key]
    return res


def generate_multilift_sequences(genomes, genome_files, sequences):
    """
    Generates multilift_sequences to be used when generating alignments.

    Args:
    - genomes: List of genome names
    - genome_files: List of genome files
    - sequences: List of sequences from user, with group assigned to each sequence [(genome, fname, seqid, group), ...]

    Returns:
    - multilift_sequences: Dictionary mapping from (genome, fname, seqid, group): SeqRecord
    """
    sequences = parse_sequence_format(sequences)

    multilift_sequences = {}
    for i, genome_file in enumerate(genome_files):
        ftype, _ = sniff_filetype(genome_file.name)
        try:
            with StringIO(genome_file.read().decode("utf-8")) as F:
                for seq in SeqIO.parse(F, ftype):
                    group = None
                    key = (genomes[i], genome_file.name, seq.id)
                    if key in sequences:
                        group = sequences[key]
                    multilift_sequences.update(
                        {(genomes[i], genome_file.name, seq.id, group): seq}
                    )
        except:
            print(f"Error reading sequences: {genome_file.name}", 3)

    return multilift_sequences


def multilift(
    uiobj_download_format,
    uiobj_uploader_alignment,
    uploaded_files,
    multilift_seq_groups,
    multilift_sequences,
    multilift_genomes,
    uiobj_aligner,
) -> None:
    multilift_download = BytesIO()
    maf_alignments = []
    igv_resources = []
    igv_genomes = []
    L = Lifter()

    with (
        tarfile.open(fileobj=multilift_download, mode="w:gz")
        if uiobj_download_format == ".tar.gz"
        else zipfile.ZipFile(multilift_download, "w")
    ) as Arc:
        if not uiobj_uploader_alignment:
            for i, seq_group in enumerate(multilift_seq_groups):
                with StringIO() as F:
                    for (
                        genome,
                        fname,
                        seqid,
                        group,
                    ), seq in multilift_sequences.items():
                        if group == seq_group:
                            SeqIO.write(
                                SeqRecord(
                                    id=f"{genome} {seqid}", description="", seq=seq.seq
                                ),
                                F,
                                "fasta",
                            )
                    returncode, result = align(F, uiobj_aligner)
                if returncode:
                    return
                aln = AlignIO.read(result, "fasta")
                L.add_alignment(aln, seq_group)

                # Write alignment as fasta
                add_to_archive(
                    Arc,
                    BytesIO(bytes(result.getvalue(), "utf-8")),
                    f"alignment/{seq_group}.fa",
                    uiobj_download_format,
                )

                # Store consensus
                cons = Seq(generate_consensus(aln))
                igv_genomes.append(SeqRecord(id=seq_group, description="", seq=cons))

                # Add to maf
                maf_alignments.append(
                    MultipleSeqAlignment(
                        [
                            SeqRecord(
                                id=f"multilift.{seq_group}", description="", seq=cons
                            )
                        ]
                        + [
                            SeqRecord(
                                id=f"{s.id}.{seq_group}", description="", seq=s.seq
                            )
                            for s in aln
                        ]
                    )
                )

        # Write consensus sequences
        with StringIO() as F:
            SeqIO.write(igv_genomes, F, "fasta")
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), "utf-8")),
                f"genome/multilift.fa",
                uiobj_download_format,
            )
        del igv_genomes

        # Write maf alignment
        with StringIO() as F:
            AlignIO.write(maf_alignments, F, "maf")
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), "utf-8")),
                f"alignment/multilift.maf",
                uiobj_download_format,
            )
            igv_resources.append(f"alignment/multilift.maf")

        # Compute maf coordinates
        for idx in range(len(maf_alignments)):
            aln = maf_alignments[idx]
            coords = [
                SeqRecord(
                    id=aln[0].id,
                    description="",
                    seq=Seq("." * aln.get_alignment_length()),
                )
            ]
            for s in aln[1:]:
                seq = list(s.seq)
                i = 0
                for j, nt in enumerate(s.seq):
                    if nt != "-":
                        seq[j] = "|" if i % 10 == 9 else "."
                        i += 1
                i = 1
                for j, nt in enumerate(seq[:]):
                    if nt == "|":
                        coord_str = str(i * 10)
                        if "-" not in seq[j : j + len(coord_str)] and j + len(
                            coord_str
                        ) < len(seq):
                            seq[j : j + len(coord_str)] = list(coord_str)
                        i += 1
                coords.append(SeqRecord(id=s.id, description="", seq=Seq("".join(seq))))
            maf_alignments[idx] = MultipleSeqAlignment(coords)

        # Write coordinates file
        with StringIO() as F:
            AlignIO.write(maf_alignments, F, "maf")
            add_to_archive(
                Arc,
                BytesIO(bytes(F.getvalue(), "utf-8")),
                f"alignment/multilift.coords.maf",
                uiobj_download_format,
            )
        igv_resources.append(f"alignment/multilift.coords.maf")
        del maf_alignments

        # Liftover data files, create annotations

        for i, genome in enumerate(multilift_genomes):
            if len(multilift_genomes) != len(uploaded_files):
                continue
            file = uploaded_files[i]
            ftype, application = sniff_filetype(file.name)
            if "data" in application:
                try:
                    new_ext, lift_file = liftover(
                        StringIO(file.read().decode("utf-8")), ftype, L, genome
                    )
                except Exception as e:
                    raise ValueError(
                        f"Error lifting over {file.name} for {genome}: {e}"
                    )
            elif "annotation" in application:
                try:
                    new_ext, lift_file = annotate(
                        StringIO(file.getvalue().decode("utf-8")),
                        L,
                        {
                            k[2]: v
                            for k, v in multilift_sequences.items()
                            if k[0] == genome
                        },
                        genome,
                    )
                except Exception as e:
                    raise ValueError(
                        f"Error applying {file.name} annotations to {genome}: {e}"
                    )
            else:
                continue
            add_to_archive(
                Arc,
                BytesIO(bytes(lift_file.getvalue(), "utf-8")),
                f"liftover/{genome}/{file.name}{new_ext}",
                uiobj_download_format,
            )
            igv_resources.append(f"liftover/{genome}/{file.name}{new_ext}")

    multilift_download.seek(0)
    return multilift_download
