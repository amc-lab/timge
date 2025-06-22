from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
import os
import mimetypes
from django.conf import settings
import numpy as np
from Bio import SeqIO
import math
import json
import math
import numpy as np
import os
from django.http import JsonResponse
from timge.utils.heatmap import generate_contact_map, generate_contact_map_bedpe
import gzip
import io
from io import BytesIO
from django.http import HttpResponse, StreamingHttpResponse
import tempfile
from pathlib import Path
import tempfile, zipfile
from timge.utils.diffStructure.main import run_differential
import pandas as pd

TRACK_ROOT_DIR = settings.TRACK_ROOT_DIR

# @csrf_exempt
# @require_http_methods(["POST"])
# def upload_tracks(request):
#     """
#     Handles the upload of track files to the server.
#     Args:
#     - request: The HTTP request containing the track files.
#     Returns:
#     - JsonResponse: A response containing the status of the upload.
#     """
#     if request.method == "POST":
#         track_files = request.FILES.getlist("track_files")
#         uuid = request.POST.get("uuid")

#         # make directory for the uuid
#         print("Track root dir:", TRACK_ROOT_DIR)
#         directory = os.path.join(TRACK_ROOT_DIR, uuid)
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         # save the files to the directory
#         for track_file in track_files:
#             print(f"Saving file: {track_file.name} to {directory}")
#             with open(os.path.join(directory, track_file.name), "wb+") as destination:
#                 for chunk in track_file.chunks():
#                     destination.write(chunk)
#         # return the directory path
#         return JsonResponse({"status": "success", "directory": directory})
#     return JsonResponse({"status": "error", "message": "Invalid request method."})

import os
import gzip
import tarfile
from io import BytesIO
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from django.http import JsonResponse
import zipfile


import os
import gzip
import tarfile
import zipfile
import shutil
import subprocess
from io import BytesIO
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods


@csrf_exempt
@require_http_methods(["POST"])
def upload_tracks(request):
    """
    Handles the upload of track files to the server.
    Extracts .gz and .tar.gz files if detected, otherwise saves as-is.
    If a file ends in .fasta or .fa, runs samtools faidx on it.
    """
    if request.method == "POST":
        track_files = request.FILES.getlist("track_files")
        uuid = request.POST.get("uuid")
        path = request.POST.get("path", [])

        if path and isinstance(path, list):
            path = "/".join(path)
        else:
            path = ""

        # make directory for the uuid
        print("Track root dir:", TRACK_ROOT_DIR)
        directory = os.path.join(TRACK_ROOT_DIR, uuid, path)
        os.makedirs(directory, exist_ok=True)

        for track_file in track_files:
            file_name = track_file.name
            file_path = os.path.join(directory, file_name)
            file_data = track_file.read()
            file_like = BytesIO(file_data)

            try:
                if file_name.endswith((".tar.gz", ".tgz")):
                    print(f"Extracting tar.gz file: {file_name}")
                    with tarfile.open(fileobj=file_like, mode="r:gz") as tar:
                        tar.extractall(path=directory)
                        print(f"Extracted contents of {file_name} to {directory}")

                elif file_name.endswith(".gz") or file_data[:2] == b"\x1f\x8b":
                    print(f"Decompressing gzip file: {file_name}")
                    with gzip.open(file_like, "rb") as gz:
                        extracted_data = gz.read()
                        base_name = os.path.splitext(file_name)[0]
                        extracted_path = os.path.join(directory, base_name)
                        with open(extracted_path, "wb") as f_out:
                            f_out.write(extracted_data)
                        print(f"Saved decompressed file to: {extracted_path}")

                        # Check if decompressed file is FASTA
                        if base_name.endswith((".fasta", ".fa")):
                            run_samtools_faidx(extracted_path)

                elif file_name.endswith(".zip"):
                    with zipfile.ZipFile(file_like, "r") as zip_ref:
                        zip_ref.extractall(directory)
                        print(f"Extracted contents of {file_name} to {directory}")
                        # You could optionally walk the extracted contents and call samtools on any .fa/.fasta files

                else:
                    print(f"Saving regular file: {file_name}")
                    with open(file_path, "wb") as f_out:
                        f_out.write(file_data)

                    if file_name.endswith((".fasta", ".fa")):
                        run_samtools_faidx(file_path)

            except Exception as e:
                return JsonResponse(
                    {
                        "status": "error",
                        "message": f"Failed to process file {file_name}: {str(e)}",
                    }
                )

        return JsonResponse({"status": "success", "directory": directory})

    return JsonResponse({"status": "error", "message": "Invalid request method."})


def run_samtools_faidx(fasta_path):
    """
    Run samtools faidx on the given FASTA file to generate an index (.fai).
    """
    try:
        print(f"Running samtools faidx on {fasta_path}")
        result = subprocess.run(
            ["samtools", "faidx", fasta_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        print(f"samtools faidx completed: {result.stdout}")
    except subprocess.CalledProcessError as e:
        print(f"samtools faidx failed: {e.stderr}")
        raise Exception(f"samtools faidx failed for {fasta_path}: {e.stderr}")


@csrf_exempt
@require_http_methods(["GET"])
def get_tracks(request):
    """
    Handles the retrieval of track files from the server.
    Args:
    - request: The HTTP request containing the uuid.
    Returns:
    - JsonResponse: A response containing the list of track files and their contents.
    """
    uuid = request.GET.get("uuid")
    if not uuid:
        return JsonResponse({"status": "error", "message": "Missing UUID."}, status=400)

    directory = os.path.join(TRACK_ROOT_DIR, uuid)
    if not os.path.exists(directory):
        return JsonResponse({"status": "success", "track_files": []})

    track_files = []
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        if not os.path.isfile(file_path):
            continue
        try:
            with open(file_path, "r") as file:
                content = file.read()
            size = os.path.getsize(file_path)
            mime_type, _ = mimetypes.guess_type(file_path)
            track_files.append(
                {
                    "name": file_name,
                    "type": mime_type or "text/plain",
                    "size": size,
                    "content": content,
                }
            )
        except Exception as e:
            continue

    return JsonResponse({"status": "success", "track_files": track_files})


@csrf_exempt
@require_http_methods(["DELETE"])
def delete_tracks(request):
    """
    Handles the deletion of track files from the server. To use when all files in an instance should be deleted.
    Args:
    - request: The HTTP request containing the uuid.
    Returns:
    - JsonResponse: A response containing the status of the deletion.
    """
    if request.method == "DELETE":
        uuid = request.GET.get("uuid")
        directory = os.path.join(TRACK_ROOT_DIR, uuid)
        if os.path.exists(directory):
            for file_name in os.listdir(directory):
                file_path = os.path.join(directory, file_name)
                os.remove(file_path)
            os.rmdir(directory)
            return JsonResponse({"status": "success", "message": "Directory deleted."})
        return JsonResponse({"status": "error", "message": "Directory not found."})
    return JsonResponse({"status": "error", "message": "Invalid request method."})


@csrf_exempt
@require_http_methods(["DELETE"])
def delete_track(request):
    """
    Handles the deletion of a specific track file from the server.
    Args:
    - request: The HTTP request containing the uuid and track name.
    Returns:
    - JsonResponse: A response containing the status of the deletion.
    """
    if request.method == "DELETE":
        uuid = request.GET.get("uuid")
        track_name = request.GET.get("track_name")
        directory = os.path.join(TRACK_ROOT_DIR, uuid)
        if os.path.exists(directory):
            file_path = os.path.join(directory, track_name)
            if os.path.exists(file_path):
                os.remove(file_path)
                return JsonResponse({"status": "success", "message": "File deleted."})
            return JsonResponse({"status": "error", "message": "File not found."})
        return JsonResponse({"status": "error", "message": "Directory not found."})
    return JsonResponse({"status": "error", "message": "Invalid request method."})


@csrf_exempt
@require_http_methods(["POST"])
def generate_heatmap(request):
    try:
        data = json.loads(request.body)
    except json.JSONDecodeError:
        return JsonResponse({"status": "error", "message": "Invalid JSON."})

    uuid = data.get("uuid")
    reads_path = data.get("file_name")
    resolution = int(data.get("resolution", 100))
    segment_1 = data.get("segment_1")
    segment_2 = data.get("segment_2")
    genome_path = data.get("genome_path")
    normalise = data.get("normalise", False)

    directory = os.path.join(TRACK_ROOT_DIR, uuid)
    reads_path = os.path.join(directory, reads_path)
    reference_path = os.path.join(directory, genome_path)
    fai_path = os.path.join(directory, genome_path + ".fai")

    if not os.path.exists(reads_path):
        return JsonResponse({"status": "error", "message": "File not found."})
    if not os.path.exists(reference_path):
        return JsonResponse({"status": "error", "message": "Reference file not found."})
    if not os.path.exists(fai_path):
        os.system(f'samtools faidx "{reference_path}"')
        print("FAI file created.")

    if reads_path.endswith(".bedpe"):
        matrix = generate_contact_map_bedpe(
            reads_path, segment_1, segment_2, fai_path, normalise, resolution
        )
    else:
        matrix = generate_contact_map(
            reads_path, segment_1, segment_2, reference_path, normalise, resolution
        )

    return JsonResponse(
        {
            "status": "success",
            "message": "Heatmap generated.",
            "matrix": matrix.tolist() if matrix is not None else None,
        }
    )


@csrf_exempt
@require_http_methods(["GET"])
def generate_fai(request):
    """
    Generates a FASTA index (FAI) file for a given reference genome.
    Args:
    - request: The HTTP request containing the uuid and genome path.
    Returns:
    - JsonResponse: A response containing the status of the generation.
    """
    uuid = request.GET.get("uuid")
    genome_path = request.GET.get("genome_path")

    directory = os.path.join(TRACK_ROOT_DIR, uuid)
    reference_path = os.path.join(directory, genome_path)

    if not os.path.exists(directory):
        return JsonResponse({"status": "error", "message": "Directory not found."})

    if not os.path.exists(reference_path):
        return JsonResponse({"status": "error", "message": "Reference file not found."})

    if os.path.exists(reference_path + ".fai"):
        return JsonResponse(
            {"status": "success", "message": "FAI file already exists."}
        )

    try:
        os.system(f'samtools faidx "{reference_path}"')
        print("FAI file created.")
        return JsonResponse({"status": "success", "message": "FAI file created."})
    except Exception as e:
        return JsonResponse({"status": "error", "message": str(e)})


def get_segment(sequence, segment):
    """
    Extracts a segment from a sequence.
    Args:
    - sequence: The sequence to extract from.
    - segment: The segment to extract.
    Returns:
    - str: The extracted segment.
    """
    print(sequence, segment)
    for record in SeqIO.parse(sequence, "fasta"):
        print(record.id, segment)
        if record.id == segment:
            return str(record.seq)

    return None


@csrf_exempt
@require_http_methods(["GET"])
def get_segments(request):
    """
    Retrieves the segments from a reference genome file.
    Args:
    - request: The HTTP request containing the uuid and genome path.
    Returns:
    - JsonResponse: A response containing the list of segments.
    """
    uuid = request.GET.get("uuid")
    genome_path = request.GET.get("genome_path")

    if not uuid or not genome_path:
        return JsonResponse(
            {"status": "error", "message": "Missing uuid or genome_path."}
        )

    directory = os.path.join(TRACK_ROOT_DIR, uuid)
    reference_path = os.path.join(directory, genome_path)

    print("Directory:", directory)
    print("Reference path:", reference_path)
    if not os.path.exists(directory):
        return JsonResponse({"status": "error", "message": "Directory not found."})

    if not os.path.exists(reference_path):
        return JsonResponse({"status": "error", "message": "Reference file not found."})

    try:
        segments = [record.id for record in SeqIO.parse(reference_path, "fasta")]
    except Exception as e:
        return JsonResponse(
            {"status": "error", "message": f"Error parsing FASTA: {str(e)}"}
        )

    return JsonResponse({"status": "success", "segments": segments})


@csrf_exempt
@require_http_methods(["POST"])
def get_files_in_path(request):
    """
    Retrieves the immediate files and folders in a given directory.
    Args:
    - request: The HTTP request containing the uuid.
    Returns:
    - JsonResponse: A response containing the list of entries.
    """
    try:
        data = json.loads(request.body)
        uuid = data.get("uuid")
        path = data.get("path", "")
    except Exception as e:
        return JsonResponse({"status": "error", "message": "Invalid JSON"}, status=400)

    if not uuid:
        return JsonResponse({"status": "error", "message": "Missing UUID."}, status=400)

    if path is not None and path != "" and path != "undefined":
        path = "/".join(path)
        print("Path:", path)
        directory = os.path.join(TRACK_ROOT_DIR, uuid, path)
    else:
        directory = os.path.join(TRACK_ROOT_DIR, uuid)

    if not os.path.exists(directory):
        return JsonResponse({"status": "success", "entries": []})

    entries = []
    for entry_name in os.listdir(directory):
        entry_path = os.path.join(directory, entry_name)

        if os.path.isdir(entry_path):
            entries.append(
                {
                    "name": entry_name,
                    "type": "directory",
                }
            )
        elif os.path.isfile(entry_path):
            try:
                size = os.path.getsize(entry_path)
                entries.append(
                    {
                        "name": entry_name,
                        "type": "file",
                        "size": size,
                    }
                )
            except Exception as e:
                continue

    return JsonResponse({"status": "success", "entries": entries})


@csrf_exempt
@require_http_methods(["POST"])
def get_files_hierarchical(request):
    """
    Recursively retrieves all files and folders under a given directory.
    Returns a full hierarchical tree.
    """
    try:
        data = json.loads(request.body)
        uuid = data.get("uuid")
        path = data.get("path", "")
    except Exception as e:
        return JsonResponse({"status": "error", "message": "Invalid JSON"}, status=400)

    if not uuid:
        return JsonResponse({"status": "error", "message": "Missing UUID."}, status=400)

    if path and isinstance(path, list):
        base_path = os.path.join(TRACK_ROOT_DIR, uuid, *path)
    else:
        base_path = os.path.join(TRACK_ROOT_DIR, uuid)

    if not os.path.exists(base_path):
        return JsonResponse({"status": "success", "entries": []})

    def build_tree(directory):
        tree = []
        try:
            for entry in os.listdir(directory):
                entry_path = os.path.join(directory, entry)
                if os.path.isdir(entry_path):
                    tree.append(
                        {
                            "name": entry,
                            "type": "directory",
                            "children": build_tree(entry_path),
                        }
                    )
                elif os.path.isfile(entry_path):
                    try:
                        size = os.path.getsize(entry_path)
                        tree.append({"name": entry, "type": "file", "size": size})
                    except Exception:
                        continue
        except Exception:
            pass
        return tree

    entries = build_tree(base_path)
    return JsonResponse({"status": "success", "entries": entries})


@csrf_exempt
@require_http_methods(["POST"])
def folder_operation(request):
    """
    Handles folder operations such as creating or deleting folders.
    Args:
    - request: The HTTP request containing the operation type, uuid, and folder name.
    Returns:
    - JsonResponse: A response containing the status of the operation.
    """
    try:
        data = json.loads(request.body)
        uuid = data.get("uuid")
        path = data.get("path", "")
        folder_name = data.get("folder_name")
        operation = data.get("operation")
    except Exception as e:
        return JsonResponse({"status": "error", "message": "Invalid JSON"}, status=400)

    if not uuid or not folder_name or not operation:
        return JsonResponse(
            {"status": "error", "message": "Missing uuid, folder_name, or operation."},
            status=400,
        )

    if path and isinstance(path, list):
        base_path = os.path.join(TRACK_ROOT_DIR, uuid, *path)
    else:
        base_path = os.path.join(TRACK_ROOT_DIR, uuid)

    if not os.path.exists(base_path):
        return JsonResponse({"status": "error", "message": "Directory not found."})

    folder_path = os.path.join(base_path, folder_name)

    if operation == "create":
        try:
            os.makedirs(folder_path)
            return JsonResponse({"status": "success", "message": "Folder created."})
        except Exception as e:
            return JsonResponse({"status": "error", "message": str(e)})

    elif operation == "delete":
        try:
            if os.path.exists(folder_path):
                os.rmdir(folder_path)
                return JsonResponse({"status": "success", "message": "Folder deleted."})
            else:
                return JsonResponse({"status": "error", "message": "Folder not found."})
        except Exception as e:
            return JsonResponse({"status": "error", "message": str(e)})

    elif operation == "rename":
        new_folder_name = data.get("new_folder_name")
        if not new_folder_name:
            return JsonResponse(
                {"status": "error", "message": "Missing new folder name."},
                status=400,
            )
        new_folder_path = os.path.join(base_path, new_folder_name)
        try:
            os.rename(folder_path, new_folder_path)
            return JsonResponse({"status": "success", "message": "Folder renamed."})
        except Exception as e:
            return JsonResponse({"status": "error", "message": str(e)})

    return JsonResponse({"status": "error", "message": "Invalid operation."})


def stream_zip(target):
    """
    Walk `target` (file or folder), write it into an in-memory ZipFile,
    and yield it in 8k chunks.
    """
    buf = BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        if os.path.isfile(target):
            # Single file: use basename
            zf.write(target, os.path.basename(target))
        else:
            # Directory: walk and zip with paths relative to the directory itself
            for root, _, files in os.walk(target):
                for fname in files:
                    full_path = os.path.join(root, fname)
                    # This makes “folder/subfolder/file.ext” inside the zip
                    rel_path = os.path.relpath(full_path, start=target)
                    zf.write(full_path, rel_path)

    buf.seek(0)
    while True:
        chunk = buf.read(8192)
        if not chunk:
            break
        yield chunk


@csrf_exempt
@require_http_methods(["GET"])
def download_path(request):
    """
    Streams back a ZIP containing either:
      - one file (if `path` is a file) or
      - the entire directory tree (if `path` is a folder).
    """
    uuid = request.GET.get("uuid")
    relpath = request.GET.get("path", "").lstrip("/")
    if not uuid:
        return JsonResponse({"status": "error", "message": "Missing uuid"}, status=400)

    base_dir = os.path.abspath(os.path.join(TRACK_ROOT_DIR, uuid))
    target = os.path.abspath(os.path.join(base_dir, relpath))

    # Security: prevent traversal outside of the uuid folder
    if not target.startswith(base_dir + os.sep):
        return JsonResponse({"status": "error", "message": "Invalid path"}, status=400)
    if not os.path.exists(target):
        return JsonResponse(
            {"status": "error", "message": "File or directory not found"}, status=404
        )

    # Decide on a sensible download‐filename
    if os.path.isdir(target):
        download_name = os.path.basename(target) or uuid
    else:
        download_name = os.path.basename(target)

    response = StreamingHttpResponse(stream_zip(target), content_type="application/zip")
    response["Content-Disposition"] = f'attachment; filename="{download_name}.zip"'
    return response


def csv_to_bedgraphs(input_csv, output_prefix=""):
    """
    Reads a CSV file with columns ['chrom', 'start_pos', 'end_pos', 'raw_pval', 'effect_size']
    and writes a separate BEDGraph file for each chromosome. Each BEDGraph will have:
        chrom  start_pos  end_pos  effect_size

    Parameters:
    - input_csv (str): Path to the input CSV file.
    - output_prefix (str): Optional prefix or directory for output files.
                           Each file will be named "<output_prefix><chrom>.bedgraph".

    Example:
        csv_to_bedgraphs("results.csv", output_prefix="bedgraphs/")
    """
    # Load the CSV into a DataFrame
    df = pd.read_csv(input_csv)

    # Ensure required columns are present
    required_cols = {"chrom", "start_pos", "end_pos", "effect_size"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")

    # Group by chromosome and write one bedgraph per group
    for chrom, group in df.groupby("chrom"):
        # Select only the start_pos, end_pos, and effect_size columns
        bedgraph_df = group[["start_pos", "end_pos", "effect_size"]].copy()

        # Sort by start_pos (BEDGraph files should be sorted)
        bedgraph_df = bedgraph_df.sort_values("start_pos")

        # Construct output filename
        filename = f"{output_prefix}{chrom}.bedgraph"

        # Write to file (no header, tab‐separated)
        bedgraph_df.to_csv(filename, sep="\t", header=False, index=False)


@csrf_exempt
@require_http_methods(["POST"])
def diff_structure(request):
    print("Received request for differential structure analysis")
    print("Request data:", request.POST, request.FILES)

    # 1) parse condition names
    conds = []
    for key, value in request.POST.items():
        if key.startswith("condition_") and key.endswith("_name"):
            idx = int(key.split("_")[1])
            while len(conds) <= idx:
                conds.append({"name": None, "files": []})
            conds[idx]["name"] = value

    # 2) collect uploaded files
    for key, f in request.FILES.items():
        parts = key.split("_")
        idx = int(parts[1])
        if idx >= len(conds) or conds[idx]["name"] is None:
            return JsonResponse(
                {"error": f"Missing name for condition index {idx}"}, status=400
            )
        conds[idx]["files"].append(f)

    # 3) write uploads to a temporary work directory
    workdir = Path(tempfile.mkdtemp())
    cond1_paths, cond2_paths = [], []
    for idx, cond in enumerate(conds):
        target_list = cond1_paths if idx == 0 else cond2_paths
        for j, upload in enumerate(cond["files"]):
            # sanitize filename
            safe_name = "".join(c for c in cond["name"] if c.isalnum() or c in "-_")
            fn = workdir / f"{safe_name}_rep{j}.txt"
            with fn.open("wb") as out:
                for chunk in upload.chunks():
                    out.write(chunk)
            target_list.append(str(fn))

    # 4) run the HMM-based differential peak caller
    out_csv = workdir / "diff_peaks.csv"
    try:
        run_differential(
            cond1_paths,
            cond2_paths,
            str(out_csv),
            n_states=2,
            cov_type="diag",
            pvalue_cutoff=0.05,
        )
    except Exception as e:
        return JsonResponse(
            {"error": "Differential analysis failed", "details": str(e)}, status=500
        )

    # 5) generate one BEDGraph per chromosome from diff_peaks.csv
    #    Each bedgraph has columns: chrom, start_pos, end_pos, effect_size
    df = pd.read_csv(out_csv)
    required_cols = {"chrom", "start_pos", "end_pos", "adj_pval"}
    missing = required_cols - set(df.columns)
    if missing:
        return JsonResponse(
            {"error": f"diff_peaks.csv is missing required columns: {missing}"},
            status=500,
        )

    bedgraph_paths = []
    for chrom, group in df.groupby("chrom"):
        bg_df = group[["chrom", "start_pos", "end_pos", "adj_pval"]].sort_values(
            "start_pos"
        )
        bedgraph_file = workdir / f"{chrom}.bedgraph"
        bg_df.to_csv(bedgraph_file, sep="\t", header=False, index=False)
        bedgraph_paths.append(bedgraph_file)

    # 6) zip results (both diff_peaks.csv and all BEDGraphs)
    zip_path = workdir / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(out_csv, arcname="diff_peaks.csv")
        for bg in bedgraph_paths:
            zf.write(bg, arcname=bg.name)

    # 7) return the ZIP as the HTTP response
    with zip_path.open("rb") as fh:
        response = HttpResponse(fh.read(), content_type="application/zip")
        response["Content-Disposition"] = (
            'attachment; filename="diff_peaks_and_bedgraphs.zip"'
        )
        return response


@csrf_exempt
@require_http_methods(["POST"])
def diff_structure_test(request):
    """
    Test endpoint for differential structure analysis.
    This is a mock implementation that simulates the process.
    """
    conds = []
    for key, value in request.POST.items():
        if key.startswith("condition_") and key.endswith("_name"):
            idx = int(key.split("_")[1])
            while len(conds) <= idx:
                conds.append({"name": None, "files": []})
            conds[idx]["name"] = value

    # 2) collect uploaded files
    for key, f in request.FILES.items():
        parts = key.split("_")
        idx = int(parts[1])
        if idx >= len(conds) or conds[idx]["name"] is None:
            return JsonResponse(
                {"error": f"Missing name for condition index {idx}"}, status=400
            )
        conds[idx]["files"].append(f)

    # 3) write uploads to a temporary work directory
    workdir = Path(tempfile.mkdtemp())
    cond1_paths, cond2_paths = [], []
    for idx, cond in enumerate(conds):
        target_list = cond1_paths if idx == 0 else cond2_paths
        for j, upload in enumerate(cond["files"]):
            # sanitize filename
            safe_name = "".join(c for c in cond["name"] if c.isalnum() or c in "-_")
            fn = workdir / f"{safe_name}_rep{j}.txt"
            with fn.open("wb") as out:
                for chunk in upload.chunks():
                    out.write(chunk)
            target_list.append(str(fn))

    # 4) run the HMM-based differential peak caller
    out_csv = workdir / "diff_peaks.csv"
    try:
        run_differential(
            cond1_paths,
            cond2_paths,
            str(out_csv),
            n_states=2,
            cov_type="diag",
            pvalue_cutoff=0.05,
        )
    except Exception as e:
        return JsonResponse(
            {"error": "Differential analysis failed", "details": str(e)}, status=500
        )

    # 5) generate one BEDGraph per chromosome from diff_peaks.csv
    #    Each bedgraph has columns: chrom, start_pos, end_pos, effect_size
    df = pd.read_csv(out_csv)
    required_cols = {"chrom", "start_pos", "end_pos", "adj_pval"}
    missing = required_cols - set(df.columns)
    if missing:
        return JsonResponse(
            {"error": f"diff_peaks.csv is missing required columns: {missing}"},
            status=500,
        )

    bedgraph_paths = []
    for chrom, group in df.groupby("chrom"):
        bg_df = group[["chrom", "start_pos", "end_pos", "adj_pval"]].sort_values(
            "start_pos"
        )
        bedgraph_file = workdir / f"{chrom}.bedgraph"
        bg_df.to_csv(bedgraph_file, sep="\t", header=False, index=False)
        bedgraph_paths.append(bedgraph_file)

    # 6) zip results (both diff_peaks.csv and all BEDGraphs)
    zip_path = workdir / "results.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(out_csv, arcname="diff_peaks.csv")
        for bg in bedgraph_paths:
            zf.write(bg, arcname=bg.name)

    # 7) return the ZIP as the HTTP response
    with zip_path.open("rb") as fh:
        response = HttpResponse(fh.read(), content_type="application/zip")
        response["Content-Disposition"] = (
            'attachment; filename="diff_peaks_and_bedgraphs.zip"'
        )
        return response


def get_fasta_file(path: str):
    """
    Reads a FASTA file and returns the first sequence record.
    Args:
        path (str): The file path to the FASTA file.
    Returns:
        SeqIO.SeqRecord: The first sequence record from the FASTA file.
    """
    with open(path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            return record
    return None


@csrf_exempt
@require_http_methods(["POST"])
def predict_rna_folds(request):
    uuid = request.POST.get("uuid")
    fasta1_path = os.path.join(TRACK_ROOT_DIR, uuid, request.POST.get("fasta1"))
    fasta2_path = os.path.join(TRACK_ROOT_DIR, uuid, request.POST.get("fasta2"))
    segment1_coords = request.POST.get("segment1_coords")
    segment2_coords = request.POST.get("segment2_coords")

    if not fasta1_path or not fasta2_path:
        return JsonResponse(
            {"error": "Both fasta1 and fasta2 must be provided."}, status=400
        )

    fasta1 = get_fasta_file(fasta1_path)
    fasta2 = get_fasta_file(fasta2_path)

    if not fasta1 or not fasta2 or not fasta1.seq or not fasta2.seq:
        return JsonResponse({"error": "Invalid or empty FASTA sequences."}, status=400)

    try:
        # Extract subsegments if specified
        if segment1_coords:
            start1, end1 = map(int, segment1_coords.split(","))
            fasta1.seq = fasta1.seq[start1:end1]
        if segment2_coords:
            start2, end2 = map(int, segment2_coords.split(","))
            fasta2.seq = fasta2.seq[start2:end2]

        seq1 = str(fasta1.seq)
        seq2 = str(fasta2.seq)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Write seq1 and seq2 to separate files
            seq1_file = tmpdir_path / "seq1.seq"
            seq2_file = tmpdir_path / "seq2.seq"
            seq1_file.write_text(seq1 + "\n")
            seq2_file.write_text(seq2 + "\n")

            duplex_input_path = tmpdir_path / "duplex_input.txt"
            with open(duplex_input_path, "w") as f:
                f.write(seq1 + "\n")
                f.write(seq2 + "\n")
                f.write("@\n")

            # Run RNAduplex with input redirected from file
            duplex_output_path = tmpdir_path / "duplex_output.txt"
            with (
                open(duplex_input_path, "r") as infile,
                open(duplex_output_path, "w") as outfile,
            ):
                subprocess.run(
                    ["RNAduplex"],
                    stdin=infile,
                    stdout=outfile,
                    stderr=subprocess.PIPE,
                    cwd=tmpdir_path,
                    check=True,
                    text=True,
                )

            if not duplex_output_path.exists():
                return JsonResponse(
                    {"error": "RNAduplex output file not found."}, status=500
                )
            print("RNAduplex output:", duplex_output_path.read_text())

            # Prepare input for RNAplot
            duplex_output_lines = duplex_output_path.read_text().strip().splitlines()
            if not duplex_output_lines:
                return JsonResponse(
                    {"error": "RNAduplex returned no output."}, status=500
                )

            # Extract only structure and sequence line
            structure_line = duplex_output_lines[0].split()[0]
            # pad left and right with dots to match the length of seq1 and seq2
            struct_sizes = structure_line.split("&")
            structure_line = (
                "." * (len(seq1) - len(struct_sizes[0]))
                + structure_line
                + "." * (len(seq2) - len(struct_sizes[1]))
            )
            rnaplot_input_path = tmpdir_path / "rnaplot_input.txt"
            rnaplot_input_path.write_text(f"{seq1}&{seq2}\n{structure_line}\n")

            # Run RNAplot using the input file
            subprocess.run(
                ["RNAplot", "-i", str(rnaplot_input_path)],
                cwd=tmpdir_path,
                check=True,
                text=True,
            )

            zip_bytes_io = BytesIO()
            with zipfile.ZipFile(zip_bytes_io, "w", zipfile.ZIP_DEFLATED) as zipf:
                for file_path in tmpdir_path.iterdir():
                    zipf.write(file_path, arcname=f"{file_path.name}")

            zip_bytes_io.seek(0)
            response = HttpResponse(zip_bytes_io.read(), content_type="application/zip")
            response["Content-Disposition"] = (
                'attachment; filename="rnafold_outputs.zip"'
            )
            return response

    except subprocess.CalledProcessError as e:
        return JsonResponse(
            {"error": "RNAduplex or RNAplot failed", "stderr": e.stderr}, status=500
        )
    except Exception as e:
        return JsonResponse(
            {"error": "Unexpected error", "details": str(e)}, status=500
        )
