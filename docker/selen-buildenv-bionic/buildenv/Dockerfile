FROM ubuntu:18.04

RUN apt update && \
  DEBIAN_FRONTEND='noninteractive' \
  DEBCONF_NONINTERACTIVE_SEEN='true' \
  apt install --yes \
  build-essential \
  gfortran \
  gmt \
  gmt-gshhg
