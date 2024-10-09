FROM ubuntu:jammy

RUN apt-get update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt-get install --yes \
    build-essential \
    openmpi-bin \
    openmpi-common \
    libopenmpi-dev \
    git


RUN useradd \
--create-home \
citcomcu_user 

USER citcomcu_user

WORKDIR /home/citcomcu_user

RUN git clone --recurse-submodules 'https://github.com/geodynamics/citcomcu.git'; cd citcomcu; make;

WORKDIR /home/citcomcu_user/citcomcu

ENV PATH="/home/citcomcu_user/citcomcu/src:${PATH}"

