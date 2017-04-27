FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

RUN pip install coverage

# update security libraries in the base image
RUN pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

###### Velvet installation
#  Directions from https://www.ebi.ac.uk/~zerbino/velvet/
#  with details at https://github.com/dzerbino/velvet/blob/master/Columbus_manual.pdf
#  https://github.com/dzerbino/velvet.git
#  Download tarball from https://www.ebi.ac.uk/~zerbino/velvet/velvet_latest.tgz
#  Each time you want to update Velvet, just use the packaged update_velvet.sh script.

WORKDIR /kb/module
RUN \
  wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_latest.tgz && \ 
  tar -zvxf velvet_latest.tgz && \
  rm -f velvet_latest.tgz && \
  cd velvet_1.2.10 && \
  #ln -s velvet_1.2.10 velvet && \
  #./update_velvet.sh && \
  make && \
  cp velvet* /kb/deployment/bin/.  

RUN mkdir -p /data && \
  cp -R velvet_1.2.10/* /data/.

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
