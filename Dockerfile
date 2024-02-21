# this builds a container that can be used to run the SurveySimulator (python and fortran)
# This container is loaded into the canfar Science Portal for use/execution.
# Can also be used directly with docker.
FROM ubuntu:latest as base
USER root
RUN apt-get update && yes | unminimize 
RUN apt-get update && yes | apt-get install wget man man-db manpages-posix git \
    build-essential zip unzip xdg-utils less emacs nano xterm vim rsync tree gfortran


# SKAHA system settings and permissions
RUN apt-get update && yes | apt-get install sssd libnss-sss libpam-sss
COPY etc/nofiles.conf /etc/security/limits.d/
COPY etc/nsswitch.conf /etc/
## see https://bugzilla.redhat.com/show_bug.cgi?id=1773148
RUN touch /etc/sudo.conf && echo "Set disable_coredump false" > /etc/sudo.conf
# generate missing dbus uuid (issue #47)
RUN dbus-uuidgen --ensure

# setup this container for skaha launching
COPY etc/startup.sh /skaha/startup.sh
RUN chmod +x /skaha/startup.sh


# setup a the needed python environment
RUN apt-get update && yes | apt-get install python3 pip
RUN pip3 install cadctap
RUN pip3 install vos
RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install 'astropy<6'
RUN pip3 install --pre astroquery
RUN pip3 install matplotlib
RUN pip3 install f90wrap
RUN pip3 install rebound
RUN pip3 install jupyter
RUN pip3 install jupyterlab

# Build the SSim
RUN mkdir /opt/SSim
RUN mkdir /opt/SSim/fortran
COPY fortran/F95 /opt/SSim/fortran/F95
COPY python /opt/SSim/python
WORKDIR /opt/SSim/fortran/F95
RUN make clean && make Driver GIMEOBJ=ReadModelFromFile
RUN cp Driver /usr/local/bin/SSim
RUN pip3 install astroplan
RUN echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | debconf-set-selections 
RUN apt-get install -y ttf-mscorefonts-installer

FROM base as deploy
WORKDIR /opt/SSim/python
RUN pip install .

# Two build sets, deploy and test
FROM base as test

RUN mkdir -p /arc/home
RUN groupadd -g 1001 testuser
RUN useradd -u 1001 -g 1001 -s /bin/bash -d /arc/home/testuser -m testuser
RUN chown -R testuser /opt/SSim
WORKDIR /opt/SSim/python
RUN pip install -e .
USER testuser
WORKDIR /arc/home/testuser
COPY etc/ReadModelFromFile.in ./
ENTRYPOINT ["/skaha/startup.sh"]
