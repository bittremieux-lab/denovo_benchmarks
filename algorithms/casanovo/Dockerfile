FROM python:3.10

RUN pip install --no-cache-dir casanovo

# Set the working directory within the image
WORKDIR /app

# Copy algorithm-related files
COPY algorithms/casanovo /app/
COPY ./*.py /app/
COPY ./*.sh /app/

# Copy dataset files
COPY datasets /app/datasets/

# Run algorithm and convert outputs
RUN chmod +x make_predictions.sh

ENTRYPOINT ["/app/make_predictions.sh"]
CMD [""]
