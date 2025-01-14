from django.urls import path

from . import views

urlpatterns = [
    path("generate_alignment/", views.generate_alignment),
]
