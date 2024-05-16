# Use Debian as base image
FROM debian:latest
# Update package lists and install necessary packages
RUN apt-get update 
RUN apt-get install -y openssl 
RUN apt-get install -y libssl-dev 
RUN apt-get install -y libglib2.0-dev 
RUN apt-get install -y libegl1-mesa
RUN apt-get install wget -y
RUN apt-get install vim -y
RUN apt-get install -y bash 

RUN chmod -R 777 /home/