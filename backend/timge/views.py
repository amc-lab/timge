from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_http_methods
import os
import mimetypes

TRACK_ROOT_DIR = "/Users/mithun/Documents/Imperial/Year 4/FYP/uploaded_data/"


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
