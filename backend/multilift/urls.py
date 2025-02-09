from django.urls import path

from . import views

urlpatterns = [
    path("multilift_sequences/", views.multilift_sequences),
    path("temp/", views.run_multilift),
    path("karyotype/", views.karyotype),
    path("circos/", views.format_circos),
]
