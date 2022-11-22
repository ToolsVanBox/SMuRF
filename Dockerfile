FROM python:3.6.1

WORKDIR /smurf

COPY requirements.txt .
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
COPY . .

ENTRYPOINT [ "python", "SMuRF.py"]

