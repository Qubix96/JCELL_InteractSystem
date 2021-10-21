#Dockerfile for JCELL® tool
FROM ubuntu:latest

WORKDIR /app
RUN apt-get update
RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8
ENV LANG en_US.utf8
RUN apt-get update
RUN apt-get install -y sudo
RUN sudo apt-get install -y git
RUN git clone https://github.com/jkubis96/JCELL.git
RUN sudo DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common
RUN sudo apt-get update

RUN sudo apt -y install default-jdk

RUN sudo apt -y install python3.8
RUN sudo apt -y install python3-pip
RUN mkdir -p JCELL/setup
RUN	mkdir -p JCELL/projects
ENV PATH="/root/miniconda3/bin:$PATH"
ARG PATH="/root/miniconda3/bin:$PATH"
RUN apt-get update
RUN pip install kedro
WORKDIR /app/JCELL
RUN	cd bin \
	&& kedro install 
ENV JUPYTER_ALLOW_INSECURE_WRITES=true
RUN	sudo apt-get install python3-knitpy -y 



CMD bash -i docker




