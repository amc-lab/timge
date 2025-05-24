from django.http import JsonResponse
from django.http import FileResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from libraries.multilift import (
    get_multilift_sequences,
    generate_multilift_sequences,
    multilift,
)
from multilift.parser import format_genome
import json
import os
from django.conf import settings
from io import BytesIO
import zipfile
import urllib.request

TRACK_ROOT_DIR = settings.TRACK_ROOT_DIR


@csrf_exempt
@require_http_methods(["POST"])
def multilift_sequences(request):
    try:
        genome = request.POST.get("genome")
        genome_file = request.FILES.get("genome_file")

        multilift_sequences = (
            {}
            if "multilift_sequences" not in request.session
            else request.session["multilift_sequences"]
        )
        get_multilift_sequences(multilift_sequences, genome_file, genome)

        return JsonResponse(list(multilift_sequences.keys()), safe=False)
    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)


@csrf_exempt
@require_http_methods(["POST"])
def run_multilift(request):
    """
    Args:
    - download_format: .zip or .tar.gz
    - genomes: list of genomes
    - genome_files: FASTA files for each genome
    - sequences: dictionary describing sequences mapping from (genome, fname, seqid) -> group
    - groups: list of groups from UI
    - aligner: alignment method to use for alignment, e.g., mafft, MUSCLE

    Returns:
    - BytesIO with alignment files
    """

    download_format = request.POST.get("download_format")
    genomes = json.loads(request.POST.get("genomes"))
    groups = json.loads(request.POST.get("groups"))
    genome_files = request.FILES.getlist("genome_files")
    uploaded_files = request.FILES.getlist("uploaded_files")
    multilift_genomes = json.loads(request.POST.get("multilift_genomes"))
    sequences = json.loads(request.POST.get("sequences"))
    aligner = request.POST.get("aligner")
    uuid = request.POST.get("uuid")

    multilift_sequences = generate_multilift_sequences(genomes, genome_files, sequences)
    print(multilift_sequences)
    res = multilift(
        download_format,
        False,
        uploaded_files,
        groups,
        multilift_sequences,
        multilift_genomes,
        aligner,
    )
    if uuid:
        path = os.path.join(TRACK_ROOT_DIR, uuid)
        os.makedirs(path, exist_ok=True)

        # Extract zip contents instead of saving the zip file
        res.seek(0)
        with zipfile.ZipFile(BytesIO(res.read())) as zip_file:
            zip_file.extractall(path)

        return JsonResponse(
            {"status": "success", "message": "Files extracted successfully"}, status=200
        )

    response = FileResponse(res, as_attachment=True, filename="multilift.zip")
    return response


@csrf_exempt
@require_http_methods(["POST"])
def karyotype(request):
    """
    Args:
    - genome: FASTA file for genome

    Returns:
    - JSON with karyotype data
    """

    genome = request.FILES.get("genome")
    data = map_genome_to_karyotype(genome)
    return JsonResponse(data, safe=False)


# @csrf_exempt
# @require_http_methods(["POST"])
# def format_circos(request):
#     """
#     Args:
#     - track_types: list of track types
#     - data_files: list of data files

#     Returns:
#     - List of dictionaries containing the circos data
#     """

#     track_types = json.loads(request.POST.get("track_types"))
#     data_files = request.FILES.getlist("data_files")
#     data = []
#     for i in range(len(track_types)):
#         data.append(format_genome(data_files[i], track_types[i]))

#     return JsonResponse(data, safe=False)


@csrf_exempt
@require_http_methods(["POST"])
def format_circos(request):
    """
    Args:
    - track_types: list of track types
    - data_files: list of data files

    Returns:
    - List of dictionaries containing the circos data
    """

    track_types = json.loads(request.POST.get("track_types"))
    track_paths = json.loads(request.POST.get("track_paths"))

    data_files = []
    for track_path in track_paths:
        if track_path.startswith("http://") or track_path.startswith("https://"):
            with urllib.request.urlopen(track_path) as response:
                f = BytesIO(response.read())
                f.name = os.path.basename(track_path)
                data_files.append(f)
        else:
            f = open(os.path.join(track_path), "rb")
            data_files.append(f)

    data = []
    for i in range(len(track_types)):
        data.append(format_genome(data_files[i], track_types[i]))

    return JsonResponse(data, safe=False)
