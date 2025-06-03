FROM gcc:latest

RUN apt-get update && apt-get install -y --no-install-recommends \
    openmpi-bin libopenmpi-dev gdb\
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

RUN make
