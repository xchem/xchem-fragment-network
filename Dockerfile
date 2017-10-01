FROM informaticsmatters/rdkit
RUN pip install tqdm nosetests
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
