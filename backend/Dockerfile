FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

RUN apt-get update && apt-get install -y \
    curl \
    build-essential \
    libpq-dev \
    wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.526-1_amd64.deb && \
    dpkg -i mafft_7.526-1_amd64.deb && \
    rm mafft_7.526-1_amd64.deb

RUN curl -sSL https://install.python-poetry.org | python3 -

ENV PATH="/root/.local/bin:$PATH"
ENV TRACK_ROOT_DIR=$TRACK_ROOT_DIR

WORKDIR /app

COPY pyproject.toml poetry.lock /app/

RUN poetry config virtualenvs.create false && poetry install --no-root --only main

COPY . /app/

EXPOSE 8000

CMD ["python", "manage.py", "runserver", "0.0.0.0:8000"]
