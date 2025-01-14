from django.http import JsonResponse
from django.http import FileResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
from libraries.multilift import run_multilift, parse_sequence
import json


@csrf_exempt
@require_http_methods(["POST"])
def generate_alignment(request):
    try:
        genomes = json.loads(request.POST.get("genomes"))
        sequence_files = request.FILES.getlist("sequences")

        multilift_sequences = {}
        for i, sequence_file in enumerate(sequence_files):
            genome = genomes[i]
            parse_sequence(genome, sequence_file, multilift_sequences)

        aligner = json.loads(request.POST.get("aligner"))
        output_format = json.loads(request.POST.get("output_format"))
        res = run_multilift(genomes, multilift_sequences, aligner, output_format)
        response = FileResponse(res, as_attachment=True, filename="multilift.zip")
        return response

    except json.JSONDecodeError:
        return JsonResponse({"error": "Invalid JSON"}, status=400)
