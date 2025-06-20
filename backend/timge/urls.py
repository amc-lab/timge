from django.urls import path

from . import views

urlpatterns = [
    path("upload_tracks/", views.upload_tracks, name="upload_tracks"),
    path("get_tracks/", views.get_tracks, name="get_tracks"),
    path("delete_tracks/", views.delete_tracks, name="delete_tracks"),
    path("delete_track/", views.delete_track, name="delete_track"),
    path("heatmap/", views.generate_heatmap, name="generate_heatmap"),
    path("get_segments/", views.get_segments, name="get_segments"),
    path("generate_fai/", views.generate_fai, name="generate_fai"),
    path("get_files_in_path/", views.get_files_in_path, name="get_files"),
    path(
        "get_files_hierarchical/",
        views.get_files_hierarchical,
        name="get_files_hierarchical",
    ),
    path("download/", views.download_path, name="download"),
    path("diff_structure/", views.diff_structure, name="diff_structure"),
    path("predict_rna_folds/", views.predict_rna_folds, name="predict_rna_folds"),
]
