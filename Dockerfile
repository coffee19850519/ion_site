FROM --platform=linux/amd64 python:3.9.13-slim

WORKDIR /app
COPY . /app

RUN pip install -r requirements.txt -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple --extra-index-url https://download.pytorch.org/whl/cu117 \
    && pip install -r requirements_torch.txt -f https://data.pyg.org/whl/torch-1.13.1+cu117.html

# Install BLAST packages
RUN apt-get update && apt-get install -y ncbi-blast+ && rm -rf /var/lib/apt/lists/*


EXPOSE 5004

CMD ["bash", "-c", "tail -f /dev/null & uvicorn main:app --host 0.0.0.0 --port 5004"]
