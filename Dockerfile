FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

#RUN apt-get update
RUN pip install --upgrade pip

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
#  defaultL: make && \ meaning 'CATEGORIES=2' and 'MAXKMERLENGTH=31'

WORKDIR /kb/module
RUN \
  wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_latest.tgz && \ 
  tar -zvxf velvet_latest.tgz && \
  ln -s velvet_1.2.10 velvet && \
  rm -f velvet_latest.tgz && \
  cd velvet && \
  #./update_velvet.sh && \
  make 'MAXKMERLENGTH=121' && \ 
  cp velvet* /kb/deployment/bin/.  

# For the testing data that comes with the software package, may not need the copying line
#RUN mkdir /velvet_data && \
  #cp -R velvet_1.2.10/data/* /velvet_data/.

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
