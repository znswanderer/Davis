FROM continuumio/miniconda
RUN conda update conda && \
    conda install numpy && \
    conda install scipy && \
    conda install pyopengl
RUN apt-get -y update && \
    apt-get -y install libgl1-mesa-glx libgl1-mesa-dri && \
    apt-get -y install mesa-utils && \
    # Install opengl a second time from apt-get because
    # in conda ther seems to be a problem with glut.
    # (When I use both it seems to work)
    apt-get -y install python-opengl && \
    apt-get -y install gcc
RUN mkdir /usr/home && \
    cd /usr/home && \
    git clone https://github.com/tscheff/Davis.git && \
    cd Davis && \
    python setup.py build_ext -i
CMD cd /usr/home/Davis && \
    python davis.py
