# FROM python:3.11-slim

# ENV PYTHONDONTWRITEBYTECODE 1
# ENV PYTHONUNBUFFERED 1

# RUN apt-get update && apt-get install -y \
#     curl \
#     build-essential \
#     libpq-dev \
#     wget && \
#     apt-get clean && rm -rf /var/lib/apt/lists/*

FROM ubuntu:24.04

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    curl build-essential libpq-dev wget \
    libgsl27 liblapack3 liblapacke libomp5 zlib1g \
    autoconf automake make gcc perl libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb && \
    dpkg -i mafft_7.526-1_amd64.deb && \
    rm mafft_7.526-1_amd64.deb

RUN apt-get update && apt-get install -y \
    libgsl27 \
    zlib1g \
    libomp5 && \
    wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_24_04/viennarna_2.7.0-1_amd64.deb && \
    dpkg -i viennarna_2.7.0-1_amd64.deb || apt-get -f install -y && \
    rm viennarna_2.7.0-1_amd64.deb

RUN curl -sSL https://install.python-poetry.org | python3 -

RUN apt-get update && \
    apt-get install -y \
        autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev \
        libcurl4-gnutls-dev libssl-dev libncurses5-dev libdeflate-dev && \
    curl -L -O https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
    tar -xjf samtools-1.21.tar.bz2 && \
    cd samtools-1.21 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.21 samtools-1.21.tar.bz2


ENV PATH="/root/.local/bin:$PATH"
ENV TRACK_ROOT_DIR=$TRACK_ROOT_DIR

WORKDIR /app

COPY pyproject.toml poetry.lock /app/

RUN poetry config virtualenvs.create false && poetry install --no-root --only main

COPY . /app/

EXPOSE 8000

CMD ["python3", "manage.py", "runserver", "0.0.0.0:8000"]
