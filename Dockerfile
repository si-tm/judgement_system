FROM ubuntu:latest

RUN apt-get update
RUN apt-get -y upgrade

RUN apt-get install python3-pip -y

RUN mkdir /home/user
COPY requirements.txt /home/user/requirements.txt

RUN pip3 install -r /home/user/requirements.txt
# RUN pip3 install -r /home/user/requirements.txt

COPY nupack-4.0.1.8/package/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl /home/user
RUN pip3 install /home/user/nupack-4.0.1.8-cp310-cp310-manylinux_2_17_aarch64.manylinux2014_aarch64.whl

# OSにあるbashを呼び出す
CMD ["/bin/bash"]