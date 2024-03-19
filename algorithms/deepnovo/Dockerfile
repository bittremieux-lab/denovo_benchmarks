FROM python:2.7

# Install required system packages
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        python-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --upgrade https://storage.googleapis.com/tensorflow/linux/cpu/tensorflow-0.10.0-cp27-none-linux_x86_64.whl

# Install Cython and pandas
RUN pip install numpy==1.16.6 pandas==0.24.2 pyteomics==4.0
RUN pip install Cython==3.0.9

# Set the working directory within the image
WORKDIR /app

# Copy algorithm-related files
COPY algorithms/deepnovo /app/
# Copy dataset files
COPY ./datasets /app/datasets/

RUN git clone -b PNAS https://github.com/nh2tran/DeepNovo.git
# Download algorithm weights
RUN curl -LRO ftp://massive.ucsd.edu/v01/MSV000081382/peak/DeepNovo/HighResolution/knapsack.npy
# Download model checkpoint (change to use other model weights)
RUN curl -LRO ftp://massive.ucsd.edu/v01/MSV000081382/peak/DeepNovo/HighResolution/train.deepnovo.high.cross.9high_80k.exclude_mouse/checkpoint
RUN curl -LRO ftp://massive.ucsd.edu/v01/MSV000081382/peak/DeepNovo/HighResolution/train.deepnovo.high.cross.9high_80k.exclude_mouse/translate.ckpt-83700

# Run algorithm and convert outputs
RUN chmod +x make_predictions.sh

ENTRYPOINT ["/app/make_predictions.sh"]
CMD [""]
