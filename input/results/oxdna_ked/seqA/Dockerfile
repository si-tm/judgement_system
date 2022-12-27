FROM ubuntu:20.04
 ENV DEBIAN_FRONTEND=noninteractive
 
RUN apt-get update && \
   apt-get install -y build-essential cmake clang python3-dev vim
RUN apt-get install -y pip
RUN pip install -U pip && \
    pip install numpy 
RUN apt update -y
RUN apt install git -y
WORKDIR code
RUN git clone --depth 1 https://github.com/si-tm/oxDNA.git
RUN mkdir oxDNA/build
RUN cd oxDNA/build && cmake .. && make -j4
RUN cd oxDNA/ && sed -i -e '/SEQ_FILE="seqA-GA100000-0.80_final_20200904131038.dat"/c SEQ_FILE="seqA-GA100000-0.80_final_20200904131038.dat"'
RUN cd oxDNA/ && chmod 755 tmp.sh && ./tmp.sh
