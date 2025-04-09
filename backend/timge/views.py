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

TRACK_ROOT_DIR = settings.TRACK_ROOT_DIR


@csrf_exempt
@require_http_methods(["POST"])
def upload_tracks(request):
    """
    Handles the upload of track files to the server.
    Args:
    - request: The HTTP request containing the track files.
    Returns:
    - JsonResponse: A response containing the status of the upload.
    """
    if request.method == "POST":
        track_files = request.FILES.getlist("track_files")
        uuid = request.POST.get("uuid")

        # make directory for the uuid
        print("Track root dir:", TRACK_ROOT_DIR)
        directory = os.path.join(TRACK_ROOT_DIR, uuid)
        if not os.path.exists(directory):
            os.makedirs(directory)
        # save the files to the directory
        for track_file in track_files:
            print(f"Saving file: {track_file.name} to {directory}")
            with open(os.path.join(directory, track_file.name), "wb+") as destination:
                for chunk in track_file.chunks():
                    destination.write(chunk)
        # return the directory path
        return JsonResponse({"status": "success", "directory": directory})
    return JsonResponse({"status": "error", "message": "Invalid request method."})


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
