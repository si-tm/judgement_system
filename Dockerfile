FROM ubuntu:latest

RUN apt-get update
RUN apt-get -y upgrade

RUN apt-get install python3-pip -y

RUN mkdir /home/user
COPY requirements.txt /home/user/requirements.txt

RUN pip3 install -r /home/user/requirements.txt
# RUN pip3 install -r /home/user/requirements.txt

# OSにあるbashを呼び出す
CMD ["/bin/bash"]