FROM informaticsmatters/rdkit
RUN pip install tqdm nose
ADD . /usr/local/fragalysis
RUN pip install /usr/local/fragalysis
