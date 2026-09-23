FROM geodynamics/selen-buildenv-bionic:latest

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
  git

RUN useradd \
  --create-home \
  selen_user

USER selen_user

WORKDIR /home/selen_user

RUN git clone 'https://github.com/geodynamics/selen.git';

RUN cd selen/DATA; gunzip *.gz;

RUN cd selen/src; make;

ENV PATH="/home/selen_user/selen:${PATH}"